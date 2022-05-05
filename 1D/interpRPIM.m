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
% ## @deftypefn {} {@var{retval} =} interpRPIM (@var{input1}, @var{input2})
% ##
% ## @seealso{}
% ## @end deftypefn
%
% ## Author: Po Chen <pochen@fresnel>
% ## Created: 2019-10-20

function [f, alpha, C0, C] = interpRPIM (x, xc, fc, shape, ch, mp, B, rb, rbx, eta, dC)

[H, re, ~] = evaluationMatrix(x, xc, shape, ch, mp);

D = mqDerivatives(rb, rbx, shape, 1, ch);
% if mp == 1
%     D = [D, zeros(length(xc), 1)];
% elseif mp == 2
%     D = [D, ones(length(xc), 1), zeros(length(xc), 1)];
% end

f1 = [fc; zeros(mp, 1)];
alpha = B \ f1;

C0 = abs(alpha(1:end-mp) .* (D * alpha(1:end-mp)));
C0 = C0 / max(C0);
C = C0;

while(true)
    
    C0 = C;
    idx = find(C0 > eta);
    B(1:end-mp, idx) = rb(:, idx);
    D(:, idx) = sign(rbx(:, idx));
    H(:, idx) = re(:, idx);
    
    alpha = B \ f1;
    f = H * alpha;
    C = abs(alpha(1:end-mp) .* (D * alpha(1:end-mp)));
    C = C / max(C);
    
    if max(abs(C - C0)) < dC
        break;
    end
    
end


