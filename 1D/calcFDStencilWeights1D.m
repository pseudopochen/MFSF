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
## @deftypefn {} {@var{retval} =} calcFDStencilWeights1D (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2019-10-15

function [Dx, Dxx] = calcFDStencilWeights1D (neighbors, x, c)
  N = length(x);
  wt = zeros(N, length(neighbors(1,:)) + 1);
  Dx = zeros(N, N);
  Dxx = zeros(N, N);
  
  for i = 1:N
    tn = neighbors(i, find(neighbors(i,:)));
    pn = [tn, i];
    o = ones(1, length(pn));
    rx = x(pn)' * o - (x(pn)' * o)';
    r = abs(rx);
    B = mq(r, c);
    %size(B)
    Hx = mqDerivatives(abs(x(i) - x(pn)), x(i) - x(pn), c, 1);
    %size(Hx)
    Hxx = mqDerivatives(abs(x(i) - x(pn)), x(i) - x(pn), c, 2);
    wt(i,:) = Hx / B;
    Dx(i, pn) = wt(i,:);
    Dxx(i, pn) = Hxx / B;
  end
  
  Dx = sparse(Dx);
  Dxx = sparse(Dxx);

endfunction
