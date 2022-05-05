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
## @deftypefn {} {@var{retval} =} derSat (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-04

function [dm, shp0] = derSpatial (pinner, ns, shp0, dshp, minK, maxK)
  
  [neighbors, pall] = getNeighbors(pinner, [], ns);

  [pn_all, shp_all] = nearOptimalShape(neighbors, pall, pinner, shp0, minK, maxK, dshp);
  
  shp0 = shp_all(end);
  
  dm = fdStencilWeight(pinner, pall, pn_all, shp_all);
  
endfunction
