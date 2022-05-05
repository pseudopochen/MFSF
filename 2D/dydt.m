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
## @deftypefn {} {@var{retval} =} dydt (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-04

function [dy, xic_all, dwc_all] = dydt (tb, y, c, neighbors_shock, shape1, nl, useCorners, corners, dshp, minK, maxK)
  %disp(['tb = ', num2str(tb)]);
  
  dy = zeros(size(y));
   
  NP = length(y) / 4;
  pnts = [y(1:NP)'; y(NP+1:2*NP)'];
  
  [nx, ny] = normalizeAll(y(2*NP+1 : 3*NP)', y(3*NP+1 : 4*NP)');
  
  %%
  
  w = shockSpeed(pnts, [nx; ny], c);
  
  dy(1:NP) = (w .* nx)';
  dy(NP+1 : 2*NP) = (w .* ny)';
  
  %%
  
  xic_all = zeros(NP, 2);
  dwc_all = zeros(NP, 1);

  shp = shape1;
  for ii = 1:NP
    %ii
    idxnbr =  neighbors_shock(ii,:);
    pc = pnts(:, idxnbr);
    wc = w(idxnbr)';
    [xic, dwc] = derCurve(pc, wc, shp, dshp, minK, maxK);
    xic_all(ii, :) = xic(nl+1,:);
    dwc_all(ii) = dwc(nl+1);
  endfor

  if useCorners
    dwc_all(corners) = 0.0;
  endif
  
  dy(2*NP+1 : 3*NP) = -dwc_all .* xic_all(:,1);
  dy(3*NP+1 : 4*NP) = -dwc_all .* xic_all(:,2);  
  
  
endfunction
