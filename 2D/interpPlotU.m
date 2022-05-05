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
## @deftypefn {} {@var{retval} =} interpPlotU (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-17

function uplt = interpPlotU (xplt, yplt, uplt, pnts, u, y)
  [M, N] = size(xplt);
  NP = length(y) / 4;
  xpoly = [y(1:NP), y(NP+1:2*NP)];
  %ns = [y(2*NP+1:3*NP), y(3*NP+1:4*NP)];
  
  xi = reshape(xplt, M*N, 1);
  yi = reshape(yplt, M*N, 1);
  
  umask = inpoly2([xi, yi], xpoly);
  
  ui = griddata(pnts(:,1), pnts(:,2), u, xi(umask), yi(umask));
  
  uu = reshape(zeros(M,N), M*N, 1);
  uu(umask) = ui;
  
  uplt = reshape(uu, M, N);
endfunction
