function [t,G,Tc,Ppv,Pl] = riyadh_day(GHI, Tamb, daylen, load_kwh)
%RIYADH_DAY  Build one representative day from measured monthly-mean data.
%   GHI      : monthly-mean daily global horizontal irradiance (kWh/m^2/day)
%   Tamb     : monthly-mean ambient temperature (degC)
%   daylen   : day length (h)
%   load_kwh : daily load energy (kWh)
%   A clear-sky-shaped profile is scaled so its daily integral equals GHI; the
%   NOCT model sets cell temperature. Returns time (h), irradiance, cell temp,
%   array PV power (kW) and load (kW).
    n=1440; t=(0:n-1)/60;  sr=12-daylen/2; ss=12+daylen/2;
    g=zeros(1,n); day=(t>sr)&(t<ss);
    g(day)=sin(pi*(t(day)-sr)/(ss-sr)).^1.4;
    Gpk=GHI*1000/(sum(g)/60);            % peak W/m^2 giving the measured daily insolation
    G=g*Gpk;  Tc=Tamb+(45-20)/800*G;
    Ppv=arrayfun(@(gg,cc) pv_array_power(gg,cc)/1000, G, Tc);
    Pl=residential_load(t, load_kwh);
end
