%% RUN_THD_SIMULINK  Build (if needed), run and analyse the Simulink inverter+LCL model.
%  Produces the grid-current waveform and its THD by FFT, and compares against the
%  design/paper value (0.76 %). Requires MATLAB + Simulink.
mdl = 'microgrid_lcl';
if ~exist([mdl '.slx'],'file')
    build_microgrid_lcl_slx;
end
d = evalin('base','LCL_design');

fprintf('Running Simulink model %s ...\n', mdl);
simOut = sim(mdl);

% --- retrieve grid current (robust to output location) ---
try
    ig = simOut.i_grid;
catch
    ig = evalin('base','i_grid');
end
t = ig.Time(:);  y = squeeze(ig.Data);  y = y(:);

% --- analyse last 3 fundamental cycles by FFT ---
dt = median(diff(t));  f0 = d.f0;
Ns = round(3/f0/dt);
s  = y(end-Ns+1:end);  ts = t(end-Ns+1:end);
win = (0.5 - 0.5*cos(2*pi*(0:Ns-1).'/(Ns-1)));
S = fft(s.*win);  Nh = floor(Ns/2)+1;  S = S(1:Nh);
freqs = (0:Nh-1).'/(Ns*dt);
[~,k1] = min(abs(freqs-f0));
lo = max(k1-2,1);  hi = min(k1+2,Nh);
A1 = sqrt(sum(abs(S(lo:hi)).^2));
mask = true(Nh,1);  mask(1:3) = false;  mask(lo:hi) = false;
THD = 100*sqrt(sum(abs(S(mask)).^2))/A1;

% --- fundamental RMS by least squares ---
w0 = 2*pi*f0;  M = [sin(w0*ts) cos(w0*ts) ones(Ns,1)];  c = M\s;
I1rms = hypot(c(1),c(2))/sqrt(2);

fprintf('\n=== Simulink LCL inverter THD ===\n');
fprintf('  Fundamental grid current : %.2f A rms (rated %.1f A)\n', I1rms, d.In_rms);
fprintf('  Grid-current THD         : %.2f %%   (paper value 0.76 %%)\n', THD);
if THD < 5, verdict = 'PASS'; else, verdict = 'FAIL'; end
fprintf('  IEEE 519 (5%%)            : %s\n', verdict);
fprintf('  LCL resonance            : %.0f Hz\n', d.fres);

% --- plots ---
figure('Name','Simulink LCL THD','Position',[100 100 950 380]);
subplot(1,2,1);
tt = (ts-ts(1))*1000;  plot(tt, s, 'LineWidth',1.4); grid on;
xlabel('Time (ms)'); ylabel('Grid current i_{L2} (A)'); title('(a) Grid current');
subplot(1,2,2);
mag = abs(S)/A1;
semilogy(freqs/1000, mag, 'LineWidth',1.2); hold on; yline(0.05,':');
xlim([0 25]); ylim([1e-4 2]); grid on;
xlabel('Frequency (kHz)'); ylabel('Magnitude (norm. to fundamental)');
title(sprintf('(b) Spectrum  (THD = %.2f%%)', THD));
