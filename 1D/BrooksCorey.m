function [fw, dfw, Krw, Kro, lw, lo] = BrooksCorey(Sw, Swr, Sws, Muw, Muo, l)

Se = (Sw - Swr) ./ (Sws - Swr);

Krw = Se.^(2/l + 3);
Kro = (1-Se).^2 .* (1-Se.^(2/l + 1));

fw = 1./(1.0 + Kro ./ Krw .* Muw ./ Muo);

lw = Krw ./ Muw;
lo = Kro ./ Muo;

dfw = -((Muo.*Muw.*power(Sw - Swr,2).*(Sw.*(2 + l + 2.*l.*power(Se,2./l)) + l.*(2.*Swr - 2.*Swr.*power(Se,2./l) - 3.*Sws) - 2.*Sws).*(Sw - Sws).* power(Se,2./l).*(Swr - Sws))./ (l.*power(Muo.*power(Sw - Swr,3).*power(Se,2./l) - Muw.*power(Sw - Sws,2).*(Swr + Sw.*power(Se,2./l) - Swr.*power(Se,2./l) - Sws),2)));
