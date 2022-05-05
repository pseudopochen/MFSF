## Copyright (C) 2022 Po Chen
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
## @deftypefn {} {@var{retval} =} dudt (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-05

function du = dudt (t, u, dm0, dm1, t0, t1, mu, u0, NP)
  %disp(['t = ', num2str(t)]); 
  
  fp0 = dm0 * u;
  
  %fp1 = dm1(1:length(u),:) * [u; ones(NP,1)*u0];
  
  %fp = fp0 + (t-t0)/(t1-t0) * (fp1 - fp0);
  fp = fp0;
  
  [~, dfdu] = fractionalFlow(u, mu);
  
  du = -dfdu .* fp;
  
endfunction
