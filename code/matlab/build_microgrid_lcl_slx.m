function mdl = build_microgrid_lcl_slx()
%BUILD_MICROGRID_LCL_SLX  Programmatically build a runnable Simulink model of the
%   single-phase full-bridge inverter + LCL filter (CORE Simulink blocks only, no
%   Simscape). Unipolar SPWM is generated with Sign blocks (double output, so no
%   data-type settings are needed), and the LCL plant is a State-Space block.
%   Bridge voltage:  vinv = (Vdc/2)*( sign(m - carrier) - sign(-m - carrier) ).
%
%   Usage:  build_microgrid_lcl_slx;  then  run_thd_simulink
%   Requires: MATLAB + Simulink (developed for R2024b; R2016b+ expected to work).

    mdl = 'microgrid_lcl';
    if bdIsLoaded(mdl), close_system(mdl,0); end
    if exist([mdl '.slx'],'file'), delete([mdl '.slx']); end
    new_system(mdl);  open_system(mdl);

    % ---------------- design parameters ----------------
    Vg_rms=230; Vg_pk=Vg_rms*sqrt(2); f0=60; w0=2*pi*f0;
    Vdc=400; fsw=10e3; Pn=6250; In_pk=(Pn/Vg_rms)*sqrt(2);
    L1=0.87e-3; L2=0.43e-3; Cf=15.67e-6; Rd=1.43; Lg=30e-6; L2t=L2+Lg;
    Ltot=L1+L2t;
    Vinv1 = Vg_pk + 1j*w0*Ltot*In_pk;   ma=abs(Vinv1)/Vdc;  phi=angle(Vinv1);

    % LCL state-space:  x=[iL1; iL2; vcf],  u=[vinv; vg],  y=iL2
    A=[ -Rd/L1,   Rd/L1,   -1/L1;
         Rd/L2t, -Rd/L2t,   1/L2t;
         1/Cf,   -1/Cf,     0    ];
    B=[ 1/L1, 0; 0, -1/L2t; 0, 0 ];
    C=[0 1 0];  D=[0 0];

    % ---------------- blocks ----------------
    addb(mdl,'simulink/Sources/Sine Wave','Mod',[30 40 70 80], ...
        'Amplitude',num2str(ma),'Frequency',num2str(w0),'Phase',num2str(phi),'Bias','0');
    addb(mdl,'simulink/Sources/Repeating Sequence','Carrier',[30 150 70 190], ...
        'rep_seq_t',['[0 ' num2str(0.5/fsw,'%.10g') ' ' num2str(1/fsw,'%.10g') ']'],'rep_seq_y','[0 1 0]');
    addb(mdl,'simulink/Math Operations/Gain','NegMod',[110 210 150 240],'Gain','-1');

    addb(mdl,'simulink/Math Operations/Sum','Sub1',[150 45 180 85],'Inputs','+-');   % m - carrier
    addb(mdl,'simulink/Math Operations/Sum','Sub2',[200 200 230 240],'Inputs','+-'); % -m - carrier
    addb(mdl,'simulink/Math Operations/Sign','Sign1',[220 50 250 80]);
    addb(mdl,'simulink/Math Operations/Sign','Sign2',[270 205 300 235]);
    addb(mdl,'simulink/Math Operations/Sum','SumS',[330 120 360 160],'Inputs','+-'); % SA - SB
    addb(mdl,'simulink/Math Operations/Gain','Vinv',[400 125 440 155],'Gain',num2str(Vdc/2));

    addb(mdl,'simulink/Sources/Sine Wave','Grid',[400 210 440 250], ...
        'Amplitude',num2str(Vg_pk),'Frequency',num2str(w0),'Phase','0','Bias','0');

    addb(mdl,'simulink/Signal Routing/Mux','Mux',[480 130 485 240],'Inputs','2');
    addb(mdl,'simulink/Continuous/State-Space','LCL',[520 150 620 210], ...
        'A',mat2str(A),'B',mat2str(B),'C',mat2str(C),'D',mat2str(D),'X0','[0;0;0]');

    addb(mdl,'simulink/Sinks/To Workspace','Igrid',[680 150 730 180], ...
        'VariableName','i_grid','SaveFormat','Timeseries','MaxDataPoints','inf');
    addb(mdl,'simulink/Sinks/Scope','Scope',[680 210 730 250]);
    addb(mdl,'simulink/Sinks/To Workspace','Vout',[680 60 730 90], ...
        'VariableName','v_inv','SaveFormat','Timeseries','MaxDataPoints','inf');

    % ---------------- connections ----------------
    add_line(mdl,'Mod/1','Sub1/1','autorouting','on');
    add_line(mdl,'Carrier/1','Sub1/2','autorouting','on');
    add_line(mdl,'Mod/1','NegMod/1','autorouting','on');
    add_line(mdl,'NegMod/1','Sub2/1','autorouting','on');
    add_line(mdl,'Carrier/1','Sub2/2','autorouting','on');
    add_line(mdl,'Sub1/1','Sign1/1','autorouting','on');
    add_line(mdl,'Sub2/1','Sign2/1','autorouting','on');
    add_line(mdl,'Sign1/1','SumS/1','autorouting','on');
    add_line(mdl,'Sign2/1','SumS/2','autorouting','on');
    add_line(mdl,'SumS/1','Vinv/1','autorouting','on');
    add_line(mdl,'Vinv/1','Mux/1','autorouting','on');
    add_line(mdl,'Grid/1','Mux/2','autorouting','on');
    add_line(mdl,'Mux/1','LCL/1','autorouting','on');
    add_line(mdl,'LCL/1','Igrid/1','autorouting','on');
    add_line(mdl,'LCL/1','Scope/1','autorouting','on');
    add_line(mdl,'Vinv/1','Vout/1','autorouting','on');

    % ---------------- solver / config ----------------
    set_param(mdl,'SolverType','Fixed-step','Solver','ode4', ...
                  'FixedStep','5e-7','StopTime','0.1','SaveOutput','off');

    assignin('base','LCL_design',struct('Vdc',Vdc,'f0',f0,'fsw',fsw,'ma',ma, ...
        'phi',phi,'In_rms',Pn/Vg_rms,'fres',1/(2*pi)*sqrt((L1+L2)/(L1*L2*Cf))));

    save_system(mdl);
    fprintf('Built and saved %s.slx (ma=%.3f, phi=%.2f deg, fres=%.0f Hz)\n', ...
            mdl, ma, phi*180/pi, 1/(2*pi)*sqrt((L1+L2)/(L1*L2*Cf)));
end

function addb(mdl, bpath, name, pos, varargin)
%ADDB  Add a block at a position with optional name/value parameters.
    add_block(bpath, [mdl '/' name], 'Position', pos, varargin{:});
end
