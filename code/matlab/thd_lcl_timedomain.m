function out = thd_lcl_timedomain()
%THD_LCL_TIMEDOMAIN  Time-domain single-phase unipolar-SPWM inverter + LCL filter.
%   Integrates the LCL state equations under 10 kHz unipolar PWM and returns the
%   grid-current THD (windowed FFT) for the LCL filter and an equal-inductance L filter.
    Vg_rms=230; Vg_pk=Vg_rms*sqrt(2); f0=60; w0=2*pi*f0;
    Vdc=400; fsw=10e3; Pn=6250; In_pk=(Pn/Vg_rms)*sqrt(2);
    L1=0.87e-3; L2=0.43e-3; Cf=15.67e-6; Rd=1.43; Lg=30e-6; L2t=L2+Lg;
    fs=2e6; dt=1/fs; Tend=6/f0; N=round(Tend/dt); tt=(0:N-1)*dt;
    Ltot=L1+L2t;
    Vinv1=Vg_pk + 1j*w0*Ltot*In_pk;  ma=abs(Vinv1)/Vdc;  phi=angle(Vinv1);
    carrier = abs(2*(tt*fsw - floor(tt*fsw+0.5)));      % 0..1 triangle
    modw = ma*sin(w0*tt + phi);
    va = Vdc*((modw >  carrier) - 0.5);                  % unipolar leg A
    vb = Vdc*(((-modw) > carrier) - 0.5);                % unipolar leg B
    vinv = va - vb;                                      % bridge output
    vg = Vg_pk*sin(w0*tt);
    % --- LCL filter (three states: iL1, iL2, vcf) ---
    iL1=0; iL2=0; vcf=0; iL2a=zeros(1,N);
    for n=1:N-1
        inode=iL1-iL2; vnode=vcf+Rd*inode;
        iL1=iL1+dt*(vinv(n)-vnode)/L1;
        iL2=iL2+dt*(vnode-vg(n))/L2t;
        vcf=vcf+dt*(inode)/Cf;
        iL2a(n+1)=iL2;
    end
    % --- L filter (equal total inductance) for comparison ---
    iL=0; iLa=zeros(1,N);
    for n=1:N-1
        iL=iL+dt*(vinv(n)-vg(n))/Ltot;  iLa(n+1)=iL;
    end
    out.THD_LCL_pct = local_thd(iL2a, tt, f0, dt);
    out.THD_L_pct   = local_thd(iLa , tt, f0, dt);
    out.I1_rms      = fundamental_rms(iL2a, tt, f0, dt);
    out.ma=ma; out.phi_deg=phi*180/pi;
    out.fres_Hz = 1/(2*pi)*sqrt((L1+L2)/(L1*L2*Cf));
    out.pass_IEEE519 = out.THD_LCL_pct < 5;
    out.wave.t=tt; out.wave.iL2=iL2a; out.wave.iL=iLa; out.wave.vg=vg;
end

function thd = local_thd(sig, tt, f0, dt)
%LOCAL_THD  Grid-current THD via windowed FFT (fundamental cluster excluded).
    Ns=round(3/f0/dt); s=sig(end-Ns+1:end);
    win=0.5-0.5*cos(2*pi*(0:Ns-1)/(Ns-1));  s=s.*win;
    S=fft(s); Nh=floor(Ns/2)+1; S=S(1:Nh);
    freqs=(0:Nh-1)/(Ns*dt); [~,k1]=min(abs(freqs-f0));
    lo=max(k1-2,1); hi=min(k1+2,Nh);
    A1=sqrt(sum(abs(S(lo:hi)).^2));            % fundamental cluster (5 bins)
    mask=true(1,Nh); mask(1:3)=false; mask(lo:hi)=false;
    Ah=sqrt(sum(abs(S(mask)).^2));             % remaining (harmonic) content
    thd=100*Ah/A1;
end

function I1rms = fundamental_rms(sig, tt, f0, dt)
    Ns=round(3/f0/dt); s=sig(end-Ns+1:end); ts=tt(end-Ns+1:end);
    w0=2*pi*f0; A=[sin(w0*ts).' cos(w0*ts).' ones(Ns,1)];
    c=A\s.'; I1rms=hypot(c(1),c(2))/sqrt(2);
end
