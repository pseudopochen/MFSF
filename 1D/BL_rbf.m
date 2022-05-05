clear all;

T = 100; % days
L = 100; % m
N = 128;
h = L / (N - 1); 
CFL = 0.5;

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
sw(1) = 0.64;

%%

addpath ~/Downloads/mqMatlabAll/

shape = 0.5;
ch = 2;
mp = 0;

[B, r, rx] = systemMatrix(x, shape, ch, mp);
Hx = mqDerivatives(r, rx, shape, 1);
dm = Hx / B;

%%

k = CFL * h;
time = 0;
while time < T
   
   if (time + k > T)
    k = T - time;
   end
   
   [F, ~, ~, ~, ~, ~] = BrooksCorey(sw, Swr, Sws, muw, muo, lam);
   F = F * vT / phi;
   dsw = -k * (dm * F);
   sw = sw + dsw;
     
   time = time + k; 
   sw(1) = Sws;
   sw(end) = Swr;
end