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
## @deftypefn {} {@var{retval} =} fdStencilWeight (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-03

function [dm] = fdStencilWeight (pinner,pall,pn_all,shp_all)
  
  x = pall(:,1);
  y = pall(:,2);
  
  N = length(pinner);
  dm = zeros(N,N);
  
  for i = 1:N
    
    pn = pn_all(i,:);
    shp = shp_all(i);
    
    o = ones(1,length(x(pn)));
    rx = x(pn)*o - (x(pn)*o)';  
    ry = y(pn)*o - (y(pn)*o)';
    r = sqrt(rx.^2 + ry.^2);

    B = mq(r,shp);
    
    [~, vx, vy] = pressure(x(i), y(i));
    
    H = vx * mqDerivatives(sqrt((x(i) - x(pn)).^2 + (y(i) - y(pn)).^2), x(i) - x(pn), shp, 1);
    H = H + vy * mqDerivatives(sqrt((x(i) - x(pn)).^2 + (y(i) - y(pn)).^2), y(i) - y(pn), shp, 1);
    
    wt = H' / B;
    dm(i, pn) = wt; 
    
  endfor
  
  dm = sparse(dm);
  
endfunction
