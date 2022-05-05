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
## @deftypefn {} {@var{retval} =} getTi (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-07

function ti = getTi (t,y)
  
  if abs(t(end)-t(end-1)) < 1e-9
    t(end-1) = [];
    y(end-1,:) = [];
  endif
  
  [NT, L] = size(y);
  NP = L / 4;
  
  ti = [];
  vi = [];
  
  ddsum = 0.0;
  
  for it = 2:NT
    dx = y(it,1:NP) - y(it-1,1:NP);
    dy = y(it,NP+1:2*NP) - y(it-1,NP+1:2*NP);
    dd = sqrt(dx.^2 + dy.^2);
    vmax = max(dd / (t(it)-t(it-1)));
    vi = [vi; vmax];
    
    dx0 = y(it-1,1:NP-1) - y(it-1,2:NP);
    dy0 = y(it-1,NP+1:2*NP-1) - y(it-1,NP+2:2*NP);
    l0 = min(sqrt(dx0.^2 + dy0.^2));
    dtt0 = l0 / vmax;
    
    dx1 = y(it,1:NP-1) - y(it,2:NP);
    dy1 = y(it,NP+1:2*NP-1) - y(it,NP+2:2*NP);
    l1 = min(sqrt(dx1.^2 + dy1.^2));
    dtt1 = l1 / vmax;
    
    ddsum = ddsum + dd;
    
    if l0 > dd
      if l0 < ddsum
        ti = [ti; t(it-1)];
        ddsum = 0;
      else
        continue;
      endif
    else
         
      dlv = (dtt1 - dtt0) / (t(it)-t(it-1));
    
      tt = t(it-1);
      ti = [ti; tt];
    
      dtt = dtt0;

      while tt + dtt <= t(it) - dtt
        tt = tt + dtt;
        ti = [ti; tt];
      
        dtt = dtt0 + dlv * (tt - t(it-1));
      endwhile
    
    endif
  
  endfor

  ti = [ti; t(end)];
endfunction
