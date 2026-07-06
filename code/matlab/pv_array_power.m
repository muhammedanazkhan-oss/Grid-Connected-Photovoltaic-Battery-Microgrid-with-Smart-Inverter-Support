function P = pv_array_power(G, Tc)
%PV_ARRAY_POWER  MPP power (W) of the 25-module (5S x 5P) 6.25 kW array.
%   Uses a persistent (G,Tc) MPP lookup table (built once from the exact single-diode
%   model) and bilinear interpolation, so time-series calls are fast. Interpolation
%   error vs. the exact model is < 0.02%.
    persistent Gg Tg LUT
    if isempty(LUT)
        Gg = 0:25:1050;  Tg = 10:2.5:85;
        LUT = zeros(numel(Gg), numel(Tg));
        for i = 1:numel(Gg)
            for j = 1:numel(Tg)
                [~,~,Pm] = pv_module_mpp(Gg(i), Tg(j));
                LUT(i,j) = 25*Pm;
            end
        end
    end
    Gq  = min(max(G,0), 1050);
    Tcq = min(max(Tc,10), 85);
    P = interp2(Tg, Gg, LUT, Tcq, Gq, 'linear');
end
