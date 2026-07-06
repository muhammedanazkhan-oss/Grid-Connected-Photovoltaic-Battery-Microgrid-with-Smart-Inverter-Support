function [res,ser] = ems_dispatch(t, Ppv, Pload, Ebat, soc_min, soc_max, eta, Pmax, slew, soc0)
%EMS_DISPATCH  Priority-dispatch EMS integrated to cyclic (periodic) steady state.
%   Surplus PV charges the battery (up to soc_max) then exports; a deficit is met
%   by discharge (down to soc_min) then grid import. slew = [] disables ramp limit.
    if nargin<10||isempty(soc0), soc0 = 0.5*(soc_min+soc_max); end
    dt = t(2)-t(1);  n = numel(t);
    soc_start = soc0;
    for day = 1:60
        soc = soc_start;
        SOC=zeros(1,n); Pb=zeros(1,n); Pgrid=zeros(1,n);
        pv2load=zeros(1,n); pv2bat=zeros(1,n); bat2load=zeros(1,n);
        grid2load=zeros(1,n); pv2grid=zeros(1,n);
        prev_pb = 0;
        for i = 1:n
            SOC(i) = soc;
            net = Ppv(i) - Pload(i);
            headroom = max(0,(soc_max-soc))*Ebat/(eta*dt);   % kW absorbable
            avail    = max(0,(soc-soc_min))*Ebat*eta/dt;     % kW deliverable
            if net >= 0
                pch = min([net, Pmax, headroom]);
                if ~isempty(slew) && prev_pb >= 0
                    pch = min(pch, prev_pb + slew*3600*dt/1000);
                end
                pdis = 0;
                pv2load(i)=Pload(i); pv2bat(i)=pch; pv2grid(i)=net-pch; Pgrid(i)=-(net-pch);
            else
                pdis = min([-net, Pmax, avail]);
                if ~isempty(slew) && prev_pb <= 0
                    pdis = min(pdis, abs(prev_pb) + slew*3600*dt/1000);
                end
                pch = 0;
                pv2load(i)=Ppv(i); bat2load(i)=pdis; grid2load(i)=(-net)-pdis; Pgrid(i)=(-net)-pdis;
            end
            Pb(i) = pch - pdis;  prev_pb = Pb(i);
            soc = soc + (eta*pch - pdis/eta)/Ebat*dt;
            soc = min(max(soc,soc_min),soc_max);
        end
        if abs(soc - soc_start) < 1e-4, break; end
        soc_start = soc;                       % converge to periodic SOC
    end
    E = @(x) sum(x)*dt;
    res.Epv=E(Ppv); res.Eload=E(Pload); res.Epv2load=E(pv2load); res.Epv2bat=E(pv2bat);
    res.Ebat2load=E(bat2load); res.Egrid_imp=E(grid2load); res.Egrid_exp=E(pv2grid);
    res.Ppv_peak=max(Ppv); res.Pimp_peak=max(max(Pgrid,0));
    res.self_suff = 100*(1 - res.Egrid_imp/res.Eload);
    res.self_cons = 100*(res.Epv2load + res.Epv2bat)/res.Epv;
    res.load_pv=100*res.Epv2load/res.Eload; res.load_bat=100*res.Ebat2load/res.Eload;
    res.load_grid=100*res.Egrid_imp/res.Eload;
    ser.t=t; ser.Ppv=Ppv; ser.Pload=Pload; ser.Pb=Pb; ser.Pgrid=Pgrid; ser.SOC=SOC;
end
