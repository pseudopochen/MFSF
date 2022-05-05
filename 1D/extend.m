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
%## @deftypefn {} {@var{retval} =} extend (@var{input1}, @var{input2})
%##
%## @seealso{}
%## @end deftypefn
%
%## Author: Po Chen <pochen@fresnel>
%## Created: 2019-11-19

function [xe, ue] = extend (x, u, h, m, BCl, ul, BCr, ur)

xl = min(x);
xr = max(x);
N = length(u);

xe = zeros(N + 2 * m, 1);
ue = zeros(N + 2 * m, 1);
q = [1:m];

xe(m-q+1) = xl - q * h;
xe(N + m + q) = xr + q * h;
xe((m+1):(N+m)) = x(1:N);

if (BCl =='P' || BCl == 'p')
  ue(m-q+1) = u(N-q);
  ue(N+m+q) = u(q+1);
  ue((m+1):(N+m)) = u(1:N);
end

if BCl == 'D'
  ue(m-q+1) = -u(q+1) + 2 * ul;
else
  ue(m-q+1) = u(q+1);
end

if BCr == 'D'
  ue(N+m+q) = -u(N-q) + 2 * ur;
else
  ue(N+m+q) = u(N-q);
end

ue((m+1):(N+m)) = u(1:N);

end
