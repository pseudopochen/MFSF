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
%## @deftypefn {} {@var{retval} =} BL_Roe (@var{input1}, @var{input2})
%##
%## @seealso{}
%## @end deftypefn
%
%## Author: Po Chen <pochen@fresnel>
%## Created: 2019-11-21

function [numflux] = BL_Roe (sl, sr, lambda, maxSat, Swr, Sws, Muw, Muo, l, vT, phi)

[fl, ~, ~, ~, ~, ~] = BrooksCorey(sl, Swr, Sws, Muw, Muo, l);
[fr, ~, ~, ~, ~, ~] = BrooksCorey(sr, Swr, Sws, Muw, Muo, l);
fl = fl * vT / phi;
fr = fr * vT / phi;

alpha = (fl - fr) ./ (sl - sr);
numflux = (alpha >= 0).* fl + (alpha < 0) .* fr;

end
