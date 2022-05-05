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
% ## @deftypefn {} {@var{retval} =} systemMatrix (@var{input1}, @var{input2})
% ##
% ## @seealso{}
% ## @end deftypefn
% 
% ## Author: Po Chen <pochen@fresnel>
% ## Created: 2019-10-20

function [B, r, rx] = systemMatrix (xc, shape, ch, mp)
[M, N] = size(xc);
if N > M
    xc = xc';
end

o = ones(1, length(xc));
rx = xc * o - (xc * o)';
r = abs(rx);

B = mq(r, shape, ch);
if mp == 1
    P = ones(length(xc), 1);
    B = [B, P; P', zeros(mp,mp)];
end

if mp == 2
    P = [xc, ones(length(xc),1)];
    B = [B, P; P', zeros(mp, mp)];
end



