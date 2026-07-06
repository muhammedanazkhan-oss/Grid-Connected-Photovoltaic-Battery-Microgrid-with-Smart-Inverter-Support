function I = pv_module_current(V, G, Tc)
%PV_MODULE_CURRENT  Terminal current of one 250 W, 60-cell module (single-diode model).
%   I = PV_MODULE_CURRENT(V,G,Tc) returns the current (A) at terminal voltage V (V),
%   in-plane irradiance G (W/m^2) and cell temperature Tc (degC).
    q = 1.602176634e-19;  kB = 1.380649e-23;
    Ns = 60; Isc = 8.68; Voc = 37.51; a = 1.18; Rs = 0.1865; Rsh = 433.6;
    ki = 0.0045; kv = -0.117; Gstc = 1000; Tstc = 25;
    if G <= 0, I = 0; return; end
    T  = Tc + 273.15;  Vt = Ns*kB*T/q;
    Iph  = (Isc + ki*(Tc-Tstc))*G/Gstc;
    VocT = Voc + kv*(Tc-Tstc);
    I0   = (Isc + ki*(Tc-Tstc)) / (exp(VocT/(a*Vt)) - 1);
    f = @(x) Iph - I0*(exp((V + x*Rs)/(a*Vt)) - 1) - (V + x*Rs)/Rsh - x;
    B = abs(Iph) + 50;              % f is monotonic in I -> this bracket always brackets the root
    I = fzero(f, [-B, B]);
end
