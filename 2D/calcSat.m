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
## @deftypefn {} {@var{retval} =} calcDerSatForT (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-05

function [u, pnts, dm0, shp0] = calcSat (t_all, p_all, dm0, u0, u, pnts, ns, shp0, dshp, minK, maxK, mu)
  
  NT = length(t_all);
  NP = length(p_all) / NT;

  %dm0 = derSpatial(pnts, ns, shp0, dshp, minK, maxK);
  
  for it = 1:NT-1
    disp(['it = ', num2str(it), ' out of ', num2str(NT-1)]);
  
    % spatial diff operator 

    pnts = [pnts; p_all(it*NP+1:(it+1)*NP,:)];
    [dm1, shp0] = derSpatial(pnts, ns, shp0, dshp, minK, maxK);
 
    % solve for u from t0 to t1
    
    tt0 = 0;
    tt1 = t_all(it+1)-t_all(it);
    
    [t1, u1] = ode45(@(t,u)dudt(t,u,dm0, dm1, tt0, tt1, mu, u0, NP), [tt0, tt1], u);
    
    % update u and pnts
    
    u = [u1(end,:)'; u0 * ones(NP,1)];
  
    dm0 = dm1;
       
  endfor %it

endfunction
