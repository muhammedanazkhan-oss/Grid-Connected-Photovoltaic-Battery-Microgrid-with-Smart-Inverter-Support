function [t,G,Tc] = clear_sky_day(peakG, sunrise, sunset, shape_exp, Tmin, Tmax, NOCT)
%CLEAR_SKY_DAY  One-minute clear-sky irradiance and NOCT cell-temperature profiles.
    if nargin<1||isempty(peakG),     peakG=1000;   end
    if nargin<2||isempty(sunrise),   sunrise=6.5;   end
    if nargin<3||isempty(sunset),    sunset=17.5;   end
    if nargin<4||isempty(shape_exp), shape_exp=1.5; end
    if nargin<5||isempty(Tmin),      Tmin=18;       end
    if nargin<6||isempty(Tmax),      Tmax=30;       end
    if nargin<7||isempty(NOCT),      NOCT=45;       end
    n = 1440;  t = (0:n-1)/60;
    G = zeros(1,n);  day = (t>sunrise) & (t<sunset);
    G(day) = peakG * sin(pi*(t(day)-sunrise)/(sunset-sunrise)).^shape_exp;
    G = max(G,0);
    Tamb = Tmin + (Tmax-Tmin) * max(sin(pi*(t-(sunrise-1))/((sunset+2)-(sunrise-1))),0);
    Tc = Tamb + (NOCT-20)/800 * G;
end
