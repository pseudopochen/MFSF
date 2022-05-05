% ## Copyright (C) 2019 Po Chen
% ##
% ## This program is free software: you can redistribute it and/or modify it
% ## under the terms of the GNU General Public License as published by
% ## the Free Software Foundation, either version 3 of the License, or
% ## (at your option) any later version.
% ##
% ## This program is distributed in the hope that it will be useful, but
% ## WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ##
% ## You should have received a copy of the GNU General Public License
% ## along with this program.  If not, see
% ## <https://www.gnu.org/licenses/>.
%
% ## -*- texinfo -*-
% ## @deftypefn {} {@var{retval} =} evaluationMatrix (@var{input1}, @var{input2})
% ##
% ## @seealso{}
% ## @end deftypefn
%
% ## Author: Po Chen <pochen@fresnel>
% ## Created: 2019-10-20

function [H, r, rx] = evaluationMatrix (x, xc, shape, ch, mp)
[m, n] = size(xc);
if n > m
    xc = xc';
end
[m, n] = size(x);
if n > m
    x = x';
end
N = length(xc);
M = length(x);
rx = x * ones(1, N) - ones(M, 1) * xc';
r = abs(rx);
H = mq(r, shape, ch);
if mp == 1
    H = [H, ones(M, 1)];
end
if mp == 2
    H = [H, x, ones(M, 1)];
end






