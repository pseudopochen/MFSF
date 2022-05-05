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
## @deftypefn {} {@var{retval} =} interpMQ (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-11

function [YI,shp,K] = interpMQ (X, Y, XI, shp, dshp, minK, maxK)
  
  [shp, B, K] = searchShape(X, shp, dshp, minK, maxK);
  
  H = evaluationMatrixMQ(X, XI, shp);
  
  m = H / B;
  
  [~, NY] = size(Y);
  [NI, ~] = size(XI);
  YI = zeros(NI, NY);
  
  for ii = 1:NY
    YI(:,ii) = m * Y(:,ii);
  endfor

endfunction
