clear all;

T = 100; % days
L = 100; % m
N = 256;
h = L / N; 
CFL = 0.010;

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

xx = [0:(h/2):L]';
x01 = [x - h / 2, x + h / 2];
x01(1,:) = [0, h / 2];
x01(end,:) = [L - h / 2, L];

%%

addpath ~/Downloads/mqMatlabAll

shape = 10;
ch = 2;
mp = 0;
[B, rb, rbx] = systemMatrix(x, shape, ch, mp);
[H, re, rex] = evaluationMatrix(xx, x, shape, ch, mp);

iB = inv(B);

Phi = mq1D_int(x01, x, shape);
iPhi = inv(Phi);

%%

FinalTime = T;
time = 0;
tstep = 0;

while(time < FinalTime)

  k = CFL * h;
  if (time + k > FinalTime)
    k = FinalTime - time;
  end
  
  sw(1) = Sws;
  sw(end) = Swr;
  
  [fw, ~, ~, ~, ~, ~] = BrooksCorey(sw, Swr, Sws, muw, muo, lam);
  %[fw1, alpha, C0, C] = interpRPIM(xx, x, fw, shape, ch, mp, B, rb, rbx, 0.5, 0.001);
  fw1 = H * iB * fw;
  
  F = [fw1(2) - fw1(1); diff(fw1(2:2:end-1)); fw1(end) - fw1(end-1)];
  dsw = -B * iPhi * F;
  sw = sw + k * dsw;
  
  time = time + k
  tstep = tstep + 1;
end
