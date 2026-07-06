function P = residential_load(t, total_kwh, base)
%RESIDENTIAL_LOAD  AC-dominated residential demand (kW) scaled to a daily energy.
    if nargin<2||isempty(total_kwh), total_kwh=54.4; end
    if nargin<3||isempty(base),      base=1.5;       end
    shape = base*ones(size(t));
    shape = shape + 0.9*exp(-((t-7.5).^2)/(2*1.1^2));    % morning
    shape = shape + 2.1*exp(-((t-15.5).^2)/(2*2.2^2));   % afternoon AC
    shape = shape + 1.9*exp(-((t-20.5).^2)/(2*1.8^2));   % evening AC
    shape = max(shape, base*0.8);
    dt = t(2)-t(1);
    P  = shape * (total_kwh / (sum(shape)*dt));
end
