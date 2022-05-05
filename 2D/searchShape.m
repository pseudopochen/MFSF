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
## @deftypefn {} {@var{retval} =} searchShape (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-11

function [shp, B, K] = searchShape (p, shp0, dshp,minK,maxK)
  
  [L, ND] = size(p);
  o = ones(1,L);
  
  if ND == 2
    x = p(:,1);
    y = p(:,2);
    rx = x*o - (x*o)';  
    ry = y*o - (y*o)';
    r = sqrt(rx.^2 + ry.^2);
  elseif ND == 1 
    r = abs(p*o - (p*o)');
  elseif ND == 3
    x = p(:,1);
    y = p(:,2);
    z = p(:,3);
    rx = x*o - (x*o)';
    ry = y*o - (y*o)';
    rz = z*o - (z*o)';
    r = sqrt(rx.^2 + ry.^2 + rz.^2);
  endif
  
  K = 1;
  c = shp0;
  while (K<minK | K>maxK)      % find a system matrix with desired K  
    B = mq(r,c);     
    K = cond(B);
    if K<minK
        c = c - dshp;
    elseif K>maxK
        c = c + dshp;
    endif   
  endwhile
  
  shp = c;

endfunction
