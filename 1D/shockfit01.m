## Copyright (C) 2021 Po Chen
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
## @deftypefn {} {@var{retval} =} shockfit01 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Po Chen <pochen@fresnel>
## Created: 2021-10-05

clear all;

L = 100; % m
A = 100; % m^2
phi = 0.2;
QT = 1; % m^3 / day
vT = QT / A;
Swr = 0.16;
Sor = 0.2;

kx = 1e-7;
mu = 1e-2;
muw = mu;
muo = mu;
lam = 2;

SBL = 0.64;
sw0 = SBL : 0.01 : 1 - Sor; sw0 = sw0';
[fw, dfw, krw, kro, lm, lo] = BrooksCorey(sw0, Swr, 1-Sor, muw, muo, lam);

t = 100;
xSw = vT / phi * t * abs(dfw);
xx = [max(xSw); L];
sw1 = [Swr; Swr];
%plot([flipud(xSw); xx], [flipud(sw0); sw1], 'k.-');

%%

xl = 0;
xs = max(xSw);
usL = SBL;
usR = Swr;
[fL, dfw, krw, kro, lm, lo] = BrooksCorey(usL, Swr, 1-Sor, muw, muo, lam);
fL = fL * vT / phi;
fR = 0;
w = (fL - fR) / (usL - usR);

%hold on;
%plot(xs, usL, 'ro')
%plot(xs, usR, 'b+')

N = length(xSw);
x = flipud(xSw);
zx = 1./(xs - xl);
z = (x - xl) * zx;

u = flipud(sw0);
U = u / zx;

%%

addpath ~/Desktop/mqMatlabAll
shape = 0.5;
o = ones(1, N);
rz = z * o - (z * o)';
r = abs(rz);
H = zeros(N, N);
%H(1,:) = mq(r(1,:), shape);
%H(2:N-1,:) = -mqDerivatives(r(2:N-1,:), rz(2:N-1,:),shape,1);
%H(N,:) = mq(r(N,:), shape);
H = -mqDerivatives(r, rz, shape, 1);
B = mq(r, shape);
dm = H / B;

%%

function fp = dF(tt, U, t,w, z, xl, Sor, SBL, Swr, muw, muo, lam, vT, phi, dm)
  xs = w * (t+tt);
  %disp(['t, tt, xs:', num2str(t), ", " , num2str(tt), ", ", num2str(xs)]);
  zx = 1.0 / (xs - xl);
  u = U * zx;
  u(1) = 1 - Sor;
  u(end) = SBL;
  [f, ~, ~, ~, ~, ~] = BrooksCorey(u, Swr, 1-Sor, muw, muo, lam);
  F = f * vT / phi - z .* u * w;
  fp = dm * F;
endfunction

%%
tic;
tlast = 600;

while t < tlast
  dt = 100;
  if t + dt > tlast
    dt = tlast - t;
  endif
  
  [t1, U1] = ode45(@(tt,U)dF(tt,U,t,w,z,xl,Sor,SBL,Swr,muw,muo,lam,vT, phi, dm), [0, dt], U);
  
  U = U1(end,:)';
 
  t = t + dt;
  xs = w * t;
  zx = 1.0 / (xs - xl);
  U(1) = (1 - Sor) / zx;
  U(N) = SBL / zx; 
endwhile
toc

plt_SF = [(xl + z /zx)/L, U * zx];
plt_SF = [plt_SF; plt_SF(end,:)];
plt_SF(end,2) = Swr;
plt_SF = [plt_SF; 1, Swr];
save('-binary', 'dump_SF.bin','plt_SF');

plot(plt_SF(:,1),plt_SF(:,2), 'k.-','LineWidth',1)

t = tlast;
sw0 = SBL : 0.01 : 1 - Sor; sw0 = sw0';
[fw, dfw, krw, kro, lm, lo] = BrooksCorey(sw0, Swr, 1-Sor, muw, muo, lam);
xSw = vT / phi * t * abs(dfw);
xx = [max(xSw); L];
sw1 = [Swr; Swr];
%plot([flipud(xSw); xx], [flipud(sw0); sw1], 'k:','LineWidth',3);