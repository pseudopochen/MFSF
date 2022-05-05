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
## @deftypefn {} {@var{retval} =} nearOptimalShape (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-03

function [pn_all, shp_all] = nearOptimalShape (st, pall, pinner, shp0, minK, maxK, dshp)

  x = pall(:,1);
  y = pall(:,2);
  
  N = length(pinner);
  shp_all = zeros(N,1);
  pn_all = [];
  
  for i=1:N
    tn = st(i,find(st(i,:)));
    pn = [tn i];                 % include the base point of the stencil
    pn_all = [pn_all; pn];
    
    shp_all(i) = searchShape([x(pn),y(pn)], shp0, dshp, minK, maxK);    
  end
endfunction
