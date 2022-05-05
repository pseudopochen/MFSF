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
## @deftypefn {} {@var{retval} =} adaptShockPnts (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2022-01-10

function [yout, resample] = adaptShockPnts (y, isClosed, useCorners, corners, thresh)
  
  NP = length(y) / 4;
  NC = length(corners);
  
  %shp = shp0;
  
  yout = zeros(size(y));
  
  for ic = 1:NC
    ib = corners(ic);
    
    if ic < NC
      ie = corners(ic+1);
      ids = ib:ie;
    else
      ie = corners(1);
      ids = [ib:NP,ie];
    endif
    l = length(ids);
    
    [yxflag, xc, yc] = swapXY(y(ids), y(ids+NP));
    
##    [shp, B] = searchShape(xc, shp, dshp, minK, maxK);
##    alpha = B \ yc;
##    
##    S = arcLength(min(xc(1), xc(end)), max(xc(1), xc(end)), alpha, xc, shp);
    
    [S, ss] = arclength(xc, yc, 'l');
    disp(['ic = ', num2str(ic), ' S = ', num2str(S)]);
    
    s = S / (l-1);
    
    resample = false;
    %ss = zeros(l,1);
    for ix = 1:l-1 %2:l
      %ss(ix) = arcLength(min(xc(1), xc(ix)), max(xc(1), xc(ix)), alpha, xc, shp);
      if abs(ss(ix) - s) > s * thresh %abs(ss(ix)-ss(ix-1) - s) > s * thresh
        resample = true;
      endif
    endfor
    
    if !resample
      yout(ids) = y(ids);
      yout(ids+NP) = y(ids+NP);
      yout(ids+2*NP) = y(ids+2*NP);
      yout(ids+3*NP) = y(ids+3*NP);
    else
      ss = [0; cumsum(ss)];
      s1 = linspace(0,S,l); s1 = s1';
      %[xc1, shp1] = interpMQ(ss, xc, s1(2:end-1), shp, dshp, minK, maxK);
      xc1 = interp1(ss, xc, s1(2:end-1), 'linear');
      %H = evaluationMatrixMQ(xc, xc1, shp);
      %yc1 = H * alpha;
      yc1 = interp1(ss, yc, s1(2:end-1), 'linear');
      xc1 = [xc(1); xc1; xc(end)];
      yc1 = [yc(1); yc1; yc(end)];
      if !yxflag
        yout(ids) = xc1; 
        yout(ids+NP) = yc1;
      else
        yout(ids) = yc1; 
        yout(ids+NP) = xc1;
      endif
      nx = y(ids + 2*NP);
      ny = y(ids + 3*NP);
      %m = H / B;
      %nx1 = m * nx;
      %ny1 = m * ny;
      nx1 = interp1(ss, nx, s1(2:end-1), 'linear');
      ny1 = interp1(ss, ny, s1(2:end-1), 'linear');
      nx1 = [nx(1); nx1; nx(end)];
      ny1 = [ny(1); ny1; ny(end)];
      yout(ids + 2*NP) = nx1;
      yout(ids + 3*NP) = ny1;
    endif
  endfor
  
  disp(['resample = ', num2str(resample)]);
 
 
endfunction
