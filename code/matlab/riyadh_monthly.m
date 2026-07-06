function A = riyadh_monthly()
%RIYADH_MONTHLY  Twelve-month analysis for Riyadh from measured monthly-mean data.
%   Monthly-mean daily GHI and ambient temperature (K.A.CARE Renewable Resource
%   Atlas / Global Solar Atlas). EDIT the two vectors below to use exact measured
%   TMY values for another Saudi site.
    mon ={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    nd  =[31 28 31 30 31 30 31 31 30 31 30 31];
    GHI =[4.0 4.9 5.8 6.5 7.4 8.0 7.7 7.2 6.5 5.5 4.4 3.8];      % kWh/m^2/day
    Tamb=[14.5 17.5 21.5 27.5 33 36 37 36.5 33 27 20.5 15.5];    % degC
    dl  =[10.6 11.2 12.0 12.9 13.5 13.8 13.7 13.1 12.3 11.4 10.7 10.4]; % h
    load=28 + 1.95*max(Tamb-18,0);                              % AC-driven daily load (kWh)
    eta=sqrt(0.95);
    fprintf('\n[4] RIYADH 12-MONTH ANALYSIS (measured monthly resource data)\n');
    fprintf('    %-4s %5s %6s %6s %7s %7s\n','Mon','GHI','PVkWh','load','SS%','SC%');
    pv=zeros(1,12); ld=zeros(1,12); imp=zeros(1,12); ex=zeros(1,12); ss=zeros(1,12); sc=zeros(1,12); bat=zeros(1,12);
    for k=1:12
        [t,~,~,Ppv,Pl]=riyadh_day(GHI(k),Tamb(k),dl(k),load(k));
        r=ems_dispatch(t,Ppv,Pl,13.5,0.40,0.80,eta,6.75,500,0.60);
        pv(k)=sum(Ppv)/60; ld(k)=sum(Pl)/60; imp(k)=r.Egrid_imp; ex(k)=r.Egrid_exp;
        ss(k)=r.self_suff; sc(k)=r.self_cons; bat(k)=r.Ebat2load;
        fprintf('    %-4s %5.1f %6.1f %6.1f %7.1f %7.1f\n',mon{k},GHI(k),pv(k),ld(k),ss(k),sc(k));
    end
    A.pv=sum(pv.*nd); A.ld=sum(ld.*nd); A.imp=sum(imp.*nd); A.ex=sum(ex.*nd); A.bat=sum(bat.*nd);
    A.ss=100*(1-A.imp/A.ld); A.sc=100*(A.pv-A.ex)/A.pv; A.ghi_yr=sum(GHI.*nd);
    fprintf('    ANNUAL: PV=%.0f kWh, load=%.0f kWh, import=%.0f, export=%.0f\n',A.pv,A.ld,A.imp,A.ex);
    fprintf('    Annual self-sufficiency=%.1f%%  PV self-consumption=%.1f%%  GHI=%.0f kWh/m2/yr\n',A.ss,A.sc,A.ghi_yr);
end
