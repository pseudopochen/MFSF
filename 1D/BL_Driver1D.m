clear all;

T = 600; % days
L = 100; % m
N = 1024;
h = L / N; 
CFL = 0.90;

A = 100; % m^2
phi = 0.2;
QT = 1; % m^3 / day
Swr = 0.16;
Sor = 0.2;
Sws = 1 - Sor;

kx = 1e-7; 
mu = 1e-2;
muw = mu;
muo = mu;
lam = 2;

vT = QT / A;
x = [0:h:L]';
sw = zeros(size(x)) + Swr; % ic

tic;
[sw] = BL_time_int('LF', x, sw, h, CFL, T, Swr, Sws, muw, muo, lam, vT, phi);
toc

plt_lf = [x/L, sw];

save('-binary','dump_LF.bin', 'plt_lf');