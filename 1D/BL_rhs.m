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
%## @deftypefn {} {@var{retval} =} BL_rhs (@var{input1}, @var{input2})
%##
%## @seealso{}
%## @end deftypefn
%
%## Author: Po Chen <pochen@fresnel>
%## Created: 2019-11-19

function [dsw] = BL_rhs (fluxType, x, sw, h, k, maxSat, Swr, Sws, Muw, Muo, l, vT, phi)

N = length(x);
[xe, se] = extend(x, sw, h, 1, 'D', Sws, 'D', Swr);

if (fluxType(1:2) == 'LF')
  dsw = -(BL_LF(se(2:N+1), se(3:N+2), 0, maxSat, Swr, Sws, Muw, Muo, l, vT, phi) - BL_LF(se(1:N), se(2:N+1), 0, maxSat, Swr, Sws, Muw, Muo, l, vT, phi)) / h;
elseif (fluxType(1:2) == 'LW')
  dsw = -(BL_LW(se(2:N+1), se(3:N+2), k / h, maxSat, Swr, Sws, Muw, Muo, l, vT, phi) - BL_LW(se(1:N), se(2:N+1), k / h, maxSat, Swr, Sws, Muw, Muo, l, vT, phi)) / h;
elseif (fluxType(1:3) == 'Roe')
  dsw = -(BL_Roe(se(2:N+1), se(3:N+2), 0, maxSat, Swr, Sws, Muw, Muo, l, vT, phi) - BL_Roe(se(1:N), se(2:N+1), 0, maxSat, Swr, Sws, Muw, Muo, l, vT, phi)) / h;
else
  disp('wrong flux type, must be LF, LW or Roe')
end

end
