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
## @deftypefn {} {@var{retval} =} getTi_1 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-10

function ti = getTi_1 (t, y, ndec)
  if abs(t(end)-t(end-1)) < 1e-9
    t(end-1) = [];
    y(end-1,:) = [];
  endif
  
  [NT, L] = size(y);
  NP = L / 4;
  
  ti = [];
  
  vv = zeros(NT-1,1);
  for it = 2:NT
    dx = y(it,1:NP) - y(it-1,1:NP);
    dy = y(it,NP+1:2*NP) - y(it-1,NP+1:2*NP);
    dd = sqrt(dx.^2 + dy.^2);
    vmax = max(dd / (t(it)-t(it-1)));
    vv(it-1) = vmax;
  endfor
  
  ll = zeros(NT,1);
  for it = 1:NT
    dx0 = y(it,1:NP-1) - y(it,2:NP);
    dy0 = y(it,NP+1:2*NP-1) - y(it,NP+2:2*NP);
    l0 = min(sqrt(dx0.^2 + dy0.^2));
    ll(it) = l0 * ndec;
  endfor
  
  tt = t(1);
  ti = [ti; tt];
  while tt < t(end)
    li = interp1(t, ll, tt);
    vi = interp1(t(1:end-1), vv, tt, 'extrap');
    dtt = li / vi;
    if tt + dtt > t(end)
      dtt1 = t(end) - tt;
      dtt0 = dtt;
      ddtt = dtt0 - dtt1;
      dtt = dtt1;
    endif
    tt = tt + dtt;
    ti = [ti; tt];
  endwhile
  
  if ddtt > dtt0 /2 && ddtt < 0.75 * dtt0
    dtc = ddtt / (length(ti)-2);
    dti = diff(ti);
    dti(1:end-1) = dti(1:end-1) - dtc;
    ticr = cumsum(dti(1:end-1));
    ti = [ti(1); ticr; ti(end)];
  endif
  if ddtt > dtt0 /2 && ddtt >= 0.75 * dtt0
    dtc = (ti(end) - ti(end-1)) / (length(ti)-2);
    dti = diff(ti);
    dti(1:end-1) = dti(1:end-1) + dtc;
    ticr = cumsum(dti(1:end-1));
    ti = [ti(1); ticr];
  endif

endfunction
