%## Copyright (C) 2019 Po Chen
%## 
%## This program is free software: you can redistribute it and/or modify it
%## under the terms of the GNU General Public License as published by
%## the Free Software Foundation, either version 3 of the License, or
%## (at your option) any later version.
%## 
%## This program is distributed in the hope that it will be useful, but
%## WITHOUT ANY WARRANTY; without even the implied warranty of
%## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%## GNU General Public License for more details.
%## 
%## You should have received a copy of the GNU General Public License
%## along with this program.  If not, see
%## <https://www.gnu.org/licenses/>.
%
%## -*- texinfo -*- 
%## @deftypefn {} {@var{retval} =} mq1D_int (@var{input1}, @var{input2})
%##
%## @seealso{}
%## @end deftypefn
%
%## Author: Po Chen <pochen@fresnel>
%## Created: 2019-11-22

function [phi] = mq1D_int (x01, xc, c)
  [n, ~] = size(x01);
  phi = zeros(n, length(xc));
  for ii = 1:n
    x0 = x01(ii, 1);
    x1 = x01(ii, 2);
    phi0 = 0.5 * (x0 - xc) .* sqrt(1 + c.^2 .* (x0 - xc).^2) + 0.5 * asinh(c .* (x0 - xc)) ./ c;
    phi1 = 0.5 * (x1 - xc) .* sqrt(1 + c.^2 .* (x1 - xc).^2) + 0.5 * asinh(c .* (x1 - xc)) ./ c;
    phi(ii,:) = phi1 - phi0;
  end
end
