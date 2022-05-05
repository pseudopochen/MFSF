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
## @deftypefn {} {@var{retval} =} generateFDStencils1D (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2019-10-15

function neighbors = generateFDStencils1D (x, ns)
  n = length(x);
  neighbors = zeros(n, ns);
  
  for i = 1:ns
    neighbors(i, :) = 1:ns; 
  end
  for i = ns + 1:n
    neighbors(i, :) = i - ns + 1 : i;
  end
  
  for i = 1:n
    if i == 1
      x0 = x(i);
    else
      x0 = x(i);
    end
    
    [r, ix] = sort(abs(x-x0));
    neighbors(i, 1:ns) = ix(2:ns+1);
  endfor

endfunction
