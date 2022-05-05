## Copyright (C) 2019 Po Chen
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
## @deftypefn {} {@var{retval} =} getUpstreamPoints (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2019-10-22

function [xi, icount] = getUpstreamPoints (v, x, dx0, B, rb, rbx, shp, ch, mp, eta, dC, dt, err)
dx = dx0;
icount = 0;  
while(true) %for ii = 1:500
v1 = interpV(v, x, dx, B, rb, rbx, shp, ch, mp, eta, dC);
dx1 = dt * v1;
%max(abs(dx1 - dx))
if max(abs(dx1 - dx)) < err
  xi = x - dx1;
  break;
end
dx = dx1;
icount = icount + 1;
end

  


