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
%## @deftypefn {} {@var{retval} =} BL_time_int (@var{input1}, @var{input2})
%##
%## @seealso{}
%## @end deftypefn
%
%## Author: Po Chen <pochen@fresnel>
%## Created: 2019-11-19

function [sw] = BL_time_int (fluxType, x, sw, h, CFL, FinalTime, Swr, Sws, Muw, Muo, l, vT, phi)
time = 0;
tstep = 0;

if (fluxType(1:2) == 'LF')
  sw(1) = 0.64;
end

while(time < FinalTime)

  maxSat = 0.0;
  k = CFL * h;

  if (fluxType(1:2) == 'LF')
    [~, dfw, ~, ~, ~, ~] = BrooksCorey(sw, Swr, Sws, Muw, Muo, l);
    maxSat = max(abs(dfw * vT / phi));
    k = k / maxSat;
  end
  
  if (time + k > FinalTime)
    k = FinalTime - time;
  end
  
  sw = sw + k * BL_rhs(fluxType, x, sw, h, k, maxSat, Swr, Sws, Muw, Muo, l, vT, phi);
  time = time + k;
  tstep = tstep + 1;
end

end
