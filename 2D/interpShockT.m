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
## @deftypefn {} {@var{retval} =} interpT (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-07

function [ti, pntsi] = interpShockT (t, y, ndec)
    
  [NT, L] = size(y);
  NP = L / 4;
  NPS = NP / ndec;
  
  ti = getTi_1(t,y,ndec);
  ni = length(ti);
  %ti = linspace(t(1),t(end),ni);
  pntsi = zeros(NPS * ni, 2);
  
  jj = 1;
  for ii = 1:ndec:NP
    xi = interp1(t,y(:,ii),ti,'spline');
    pntsi(jj:NPS:end,1) = xi;
    jj = jj+1;  
  endfor
  
  jj = 1;
  for ii = 1:ndec:NP
    yi = interp1(t,y(:,ii+NP),ti,'spline');
    pntsi(jj:NPS:end,2) = yi;
    jj = jj + 1;
  endfor
  
endfunction
