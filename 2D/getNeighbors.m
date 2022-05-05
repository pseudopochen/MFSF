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
## @deftypefn {} {@var{retval} =} getNeighbors (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-04

function [neighbors, P] = getNeighbors (pinner, pbndy, ns)
  
  n = length(pinner);
  neighbors = zeros(n, ns);
  
  P = pinner;
  if !isempty(pbndy)
    P = [P; pbndy];
  endif

  for ii = 1:n
    x0 = pinner(ii,1);
    y0 = pinner(ii,2);
    
    r2 = (P(:,1)-x0).^2 + (P(:,2)-y0).^2;
    [r2,ix]=sort(r2);
    neighbors(ii,:) = ix(2:ns+1);
  endfor
  
endfunction
