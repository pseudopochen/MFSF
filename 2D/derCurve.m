## Copyright (C) 2022 Po Chen
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} interp2D (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-02

function [xic, dwc, shape] = derCurve (p, w, shp, dshp, minK, maxK)

[yxflag,xc,yc] = swapXY(p(1,:)', p(2,:)');

n = length(xc);

o = ones(1, n);
rx = xc * o - (xc*o)';
r = abs(rx);

a = min(abs(diff(xc)));
shape = shp / a;
B = mq(r, shape);
K = cond(B);
if K < minK || K > maxK
    [shape, B, ~] = searchShape([xc,yc], shape, dshp, minK, maxK);
    %disp(['K = ', num2str(K), ' shape = ', num2str(shape)]);
endif
Hx = mqDerivatives(r, rx, shape, 1);
shape = shape * a;

D = Hx / B;

dyc = D * yc;

if !isempty(w)
  dwc = D * w;
end

xic = zeros(n,2);
xinorm = sqrt(1 + dyc.^2);
for ii = 1:n
  if yxflag
    xi = [dyc(ii), 1] / xinorm(ii);
  else
    xi = [1, dyc(ii)] / xinorm(ii);
  endif
  xic(ii,:) = xi;
endfor

dwc = dwc ./ xinorm;

endfunction
