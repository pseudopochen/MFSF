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
## @deftypefn {} {@var{retval} =} swapXY (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-10

function [yxflag, xc1, yc1] = swapXY (xc, yc)
  
  yxflag = false;
  xc1 = xc;
  yc1 = yc;
  
  if max(xc) - min(xc) < (max(yc) - min(yc))
    %disp(['switch x and y axes']);
    yxflag = true;
    xc1 = yc;
    yc1 = xc;
  endif

endfunction
