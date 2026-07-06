function [Vmpp,Impp,Pmpp,Voc,Isc] = pv_module_mpp(G, Tc)
%PV_MODULE_MPP  Maximum-power point of one module by scanning the I-V curve.
    kv = -0.117; Tstc = 25; Voc0 = 37.51;
    if G <= 0, Vmpp=0; Impp=0; Pmpp=0; Voc=0; Isc=0; return; end
    VocT = Voc0 + kv*(Tc-Tstc);
    Vs = linspace(0, max(VocT,1), 400);
    Is = arrayfun(@(v) max(pv_module_current(v,G,Tc),0), Vs);
    Ps = Vs .* Is;
    [Pmpp, j] = max(Ps);
    Vmpp = Vs(j);  Impp = Is(j);
    Isc  = max(pv_module_current(0,G,Tc), 0);
    Voc  = VocT;
end
