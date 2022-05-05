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
## @deftypefn {} {@var{retval} =} evaluationMatrixMQ (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-11

function H = evaluationMatrixMQ (X, XI, shp)
  
  [N, ND] = size(X);
  [M, ~] = size(XI);
  
  oc = ones(1, N);
  oi = ones(M, 1);
  if ND == 1
    r = abs(XI * oc - oi * X');
  elseif ND == 2
    xcx = X(:,1);
    xcy = X(:,2);
    xx = XI(:,1);
    xy = XI(:,2);
    rx = xx * oc - oi * xcx';
    ry = xy * oc - oi * xcy';
    r = sqrt(rx.^2 + ry.^2);
  elseif ND == 3 
    xcx = X(:,1);
    xcy = X(:,2);
    xcz = X(:,3);
    xx = XI(:,1);
    xy = XI(:,2);
    xz = XI(:,3);
    rx = xx * oc - oi * xcx';
    ry = xy * oc - oi * xcy';
    rz = xz * oc - oi * xcz';
    r = sqrt(rx.^2 + ry.^2 + rz.^2);
  endif
  
  H = mq(r, shp);
endfunction
