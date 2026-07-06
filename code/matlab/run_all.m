%% RUN_ALL  Reproduce the grid-connected PV-battery microgrid study (solar-4407776).
%  Single-diode PV, 24 h EMS dispatch, time-domain LCL THD, seasonal/sensitivity,
%  and techno-economics. Prints a summary and (if plotting is available) saves figures.
%  Runs in MATLAB R2024b (base) or GNU Octave. No toolboxes required.
clear; clc;
fprintf('\n==============================================================\n');
fprintf(' PV-Battery Microgrid - MATLAB verification of solar-4407776\n');
fprintf('==============================================================\n');

%% 1. Single-diode PV verification
[Vm,Im,Pm] = pv_module_mpp(1000,25);
Parr_stc = 25*Pm/1000;
P30 = pv_array_power(1000,30)/1000;
P45 = pv_array_power(1000,45)/1000;
P25 = pv_array_power(1000,25)/1000;
derate = (P25-P45)/P25*100;
fprintf('\n[1] PHOTOVOLTAIC MODEL (single-diode)\n');
fprintf('    Module MPP (STC)      : %.1f W  (Vmpp=%.2f V, Impp=%.2f A)\n',Pm,Vm,Im);
fprintf('    Array MPP (STC)       : %.2f kW  (nameplate error %.2f%%)\n',Parr_stc,abs(Parr_stc*1000-6250)/6250*100);
fprintf('    Array @30C / @45C     : %.2f / %.2f kW   derating 25->45C = %.2f%%\n',P30,P45,derate);

%% 2. Reference-day 24 h dispatch (baseline & optimised)
[t,G,Tc,Ppv,Pload] = riyadh_day(5.98, 26.6, 12.1, 44.8);   % Riyadh annual-mean day
eta = sqrt(0.95);   % 0.9747 one-way (95% round-trip)
[b,bs] = ems_dispatch(t,Ppv,Pload,13.5,0.50,0.80,eta,6.75,[]  ,0.65);
[o,os] = ems_dispatch(t,Ppv,Pload,13.5,0.40,0.80,eta,6.75,500 ,0.60);
fprintf('\n[2] 24 h DISPATCH  (daily PV %.1f kWh, PV peak %.2f kW, load %.1f kWh)\n',sum(Ppv)/60,max(Ppv),sum(Pload)/60);
fprintf('    %-28s  %8s  %8s\n','metric','baseline','optimised');
fprintf('    %-28s  %8.1f  %8.1f\n','Self-sufficiency (%)',b.self_suff,o.self_suff);
fprintf('    %-28s  %8.1f  %8.1f\n','PV self-consumption (%)',b.self_cons,o.self_cons);
fprintf('    %-28s  %8.1f  %8.1f\n','Grid import (kWh)',b.Egrid_imp,o.Egrid_imp);
fprintf('    %-28s  %8.1f  %8.1f\n','Battery->load (kWh)',b.Ebat2load,o.Ebat2load);

%% 3. Grid-current THD (time-domain PWM + LCL + FFT)
thd = thd_lcl_timedomain();
fprintf('\n[3] GRID-CURRENT THD (time-domain, unipolar SPWM 10 kHz)\n');
fprintf('    Fundamental current   : %.1f A (rated 27.2 A), m=%.3f\n',thd.I1_rms,thd.ma);
fprintf('    THD  LCL / L filter   : %.2f%% / %.1f%%   (IEEE 519 5%%: %s)\n',thd.THD_LCL_pct,thd.THD_L_pct,ternary(thd.pass_IEEE519,'PASS','FAIL'));
fprintf('    LCL resonance         : %.0f Hz\n',thd.fres_Hz);

%% 4. Riyadh 12-month analysis (measured monthly resource data)
A = riyadh_monthly();

%% 5. Sensitivity (battery size, fixed reference load)
fprintf('\n[5] BATTERY-SIZE SENSITIVITY (reference 54.4 kWh load)\n');
fprintf('    %8s %10s %10s\n','Ebat kWh','SS %','SC %');
for E=[6.75 10 13.5 20 27]
    r=ems_dispatch(t,Ppv,Pload,E,0.40,0.80,eta,6.75,500,0.55);
    fprintf('    %8.2f %10.1f %10.1f\n',E,r.self_suff,r.self_cons);
end

%% 6. Techno-economics (seasonally weighted annual)
eco = technoeconomics(A);
fprintf('\n[6] TECHNO-ECONOMICS  (Saudi Riyal; SEC residential tariff; 1 USD = 3.75 SAR)\n');
fprintf('    CAPEX   : SAR %.0f  (~$%.0f)     LCOE SAR %.3f/kWh (~$%.3f)\n',eco.capex,eco.capex/3.75,eco.lcoe,eco.lcoe/3.75);
fprintf('    Payback (0.18 SAR/kWh): %.1f yr   NPV(25yr) SAR %.0f (~$%.0f)\n',eco.payback,eco.npv,eco.npv/3.75);
fprintf('    Battery: %.0f EFC/yr, calendar life ~%.0f yr\n',eco.efc,eco.batt_life);
fprintf('    Tariff sensitivity (SAR/kWh -> payback yr):\n');
for i=1:numel(eco.tariffs)
    fprintf('        %.2f SAR -> %.1f yr  (NPV SAR %.0f)\n',eco.tariffs(i),eco.paybacks(i),eco.npvs(i));
end

fprintf('
[7] Done. (Figure generation is intentionally omitted from this archive.)
');

function s = ternary(c,a,b)
    if c, s=a; else, s=b; end
end
