clear all;

addpath ~/Downloads/mqMatlabAll

T = 100; % days
dt = 0.125; % day
dx = 0.125; % m
ns = 8; 
c = 0.5;
%shape = 3;
visco = 1e-3;
err = 1e-6;


L = 100; % m
A = 100; % m^2
phi = 0.2;
QT = 1; % m^3 / day
Swr = 0.16;
Sor = 0.2;

kx = 1e-7; 
mu = 1e-2;
muw = mu;
muo = mu;
lam = 2;

%%

x = 0:dx:L; x = x';
N = length(x);

Sw = zeros(size(x)) + Swr; % ic
Sw(1) = 1 - Sor; % left bc
Sw(end) = Swr; % right bc

[~, dfw, ~, ~, lw, lo] = BrooksCorey(Sw, Swr, 1-Sor, muw, muo, lam);
ll = lw + lo;

dpdx = -QT / A / kx / lw(1);
p = x * dpdx - x(end) * dpdx;

neighbors = generateFDStencils1D(x, ns);
[Dx, Dxx] = calcFDStencilWeights1D(neighbors, x', c);
DD = Dxx * visco;

HS = phi / dt * eye(size(DD)) - DD;
%HS(1,:) = 0.0;
%HS(1,1) = 1.0;
%HS(end,:) = 0.0;
%HS(end,end) = 1.0;
%HSinv = inv(HS);

%o = ones(1,length(x));
%rx = x * o - (x * o)';
%r = abs(rx);
%B = mq(r, shape);
%Hx = mqDerivatives(r, rx, shape, 1);
%D = Hx / B;

%%

Sw0 = Sw;
Sw1 = Sw;
p0 = p;
p1 = p;

for t = dt : dt : T
  disp(['t = ', num2str(t)]);

  icount = 0;
  do %for jj = 1:1 
    Sw = Sw1;
    p = p1;

    %dSw = dt / phi * kx * Dx * diag(lw) * Dx * p1;
    %Sw1 = Sw0 + dSw;
    %Sw1(1) = 1 - Sor;
    %Sw1(end) = Swr;
    
    %f = phi / dt * Sw0 + kx * Dx * diag(lw) * Dx * p1;
    %f(1) = 1 - Sor;
    %f(end) = Swr;
    %Sw1 = HSinv * f;
    
    vT = -kx * diag(ll) * Dx * p1;
    vTdfw = vT .* dfw;
    HS1 = HS + diag(vTdfw) * Dx;
    HS1(1,:) = 0;
    HS1(1,1) = 1.0;
    HS1(end,:) = 0;
    HS1(end,end) = 1.0;
    f = phi / dt * Sw0;
    f(1) = 1 - Sor;
    f(end) = Swr;
    Sw1 = HS1 \ f;
    
    
    [~, dfw, ~, ~, lw, lo] = BrooksCorey((Sw0 + Sw1) / 2.0, Swr, 1-Sor, muw, muo, lam);
    ll = lw + lo;
    
    H = Dx * diag(ll) * Dx;
    H(1,:) = Dx(1,:);
    H(end,:) = 0.0;
    H(end,end) = 1;
    b = zeros(N, 1);
    b(1) = -QT / lw(1) / A / kx;
    p1 = H \ b;
    
    icount = icount + 1;
    %disp(['icount: ', num2str(icount), ' ', num2str(max(abs(Sw - Sw1))), ' ', num2str(max(abs(p - p1)))]);
  until (max(abs(Sw - Sw1)) < err && max(abs(p - p1)) < err)
  %end

  Sw0 = Sw1;
  p0 = p1;
end
