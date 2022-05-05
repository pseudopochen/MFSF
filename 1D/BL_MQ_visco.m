clear all;

addpath ~/Downloads/mqMatlabAll

T = 100; % days
dt = 1; % day
dx = 0.125; % m
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

Sw = zeros(size(x)) + Swr; % ic
Sw(1) = 1 - Sor; % left bc
Sw(end) = Swr; % right bc

[~, dfw, ~, ~, lw, lo] = BrooksCorey(Sw, Swr, 1-Sor, muw, muo, lam);
ll = lw + lo;

dpdx = -QT / A / kx / lw(1);
p = x * dpdx - x(end) * dpdx;

%%
N = length(x);

cMin = 1.8 / 1.5;
cMax = 2.5 / 1.5;
cMin = cMin^2;                              
cMax = cMax^2; 
c = cMin*(cMax/cMin).^((0:N-1)./(N-1)); %linspace(0.05, 50, length(x));
shape = repmat(c, [N, 1]);
%shape = 4;
%shape = 2 * ones(N,N);
%shape(:,1:10) = 3;

o = ones(1,length(x));
rx = x * o - (x * o)';
r = abs(rx);

B = mq(r, shape);
Hx = mqDerivatives(r, rx, shape, 1);
Hxx = mqDerivatives(r, rx, shape, 2);
D = Hx / B;
DD = Hxx / B;
%DD = DD * DD * DD;

visco = 2.5e-3;
DD = DD * visco;

HS = phi / dt * eye(size(DD)) - DD;
HS(1,:) = 0.0;
HS(1,1) = 1.0;
HS(end,:) = 0.0;
HS(end,end) = 1.0;
HSinv = inv(HS);

%%

Sw0 = Sw;
Sw1 = Sw;
p0 = p;
p1 = p;
%q = D * p1;

for t = dt : dt : T
  disp(['t = ', num2str(t)]);

  icount = 0;
  do 
    Sw = Sw1;
    p = p1;
    
    %%dSw = dt * kx / phi * D * diag(lw) * D * p1; dSw = dSw(2:end-1);
    %dSw = dt * kx / dx^2 / phi * (lw(2:end-1).*p1(3:end) - (lw(2:end-1) + lw(1:end-2)) .* p1(2:end-1) + lw(1:end-2).*p1(1:end-2));
    %%pause
    %Sw1(2:end-1) = Sw0(2:end-1) + dSw;
    %Sw1(1) = 1 - Sor;
    %Sw1(end) = Swr;

    %vT = -kx * diag(ll) * D * p1;
    %vTdfw = vT .* dfw;
    %HS1 = HS + diag(vTdfw) * D;
    %HS1(1,:) = 0;
    %HS1(1,1) = 1.0;
    %HS1(end,:) = 0;
    %HS1(end,end) = 1.0;
    
    f = phi / dt * Sw0 + kx * D * diag(lw) * D * p1;
    f(1) = 1 - Sor;
    f(end) = Swr;
    %Sw1 = HS1 \ f;
    Sw1 = HSinv * f;
    %Sw1 = HS \ f;
    
    %pause
    
    [~, dfw, ~, ~, lw, lo] = BrooksCorey((Sw0 + Sw1) / 2.0, Swr, 1-Sor, muw, muo, lam);
    ll = lw + lo;

    H = D * diag(ll) * D;
    H(1,:) = D(1,:);
    H(end,:) = 0.0;
    H(end,end) = 1;
    b = zeros(N, 1);
    b(1) = -QT / lw(1) / A / kx;
    p1 = H \ b;
    
    %H1 = Hx;
    %H1(1,:) = B(1,:);
    %f = zeros(N, 1);
    %f(1) = -QT / lw(1) / A / kx * ll(1);
    %q = B * (H1 \ f) ./ ll;
    
    %H1 = Hx;
    %H1(end,:) = B(end,:);
    %q(end) = 0.0;
    %p1 = B * (H1 \ q);

    %pause
    
    icount = icount +1;
    %disp(['icount: ', num2str(icount), ' ', num2str(max(abs(Sw - ...
    %Sw1))), ' ', ...
    %num2str(max(abs(p - p1)))]);

  until (max(abs(Sw - Sw1)) < err && max(abs(p - p1)) < err)

  Sw0 = Sw1;
  p0 = p1;
end