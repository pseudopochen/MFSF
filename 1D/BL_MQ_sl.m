clear all;

addpath ~/Downloads/mqMatlabAll

T = 100; % days
dt = 1; % day
dx = 2; % m
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

shp = 0.06;
mp = 0;
ch = 2;
eta = 0.5;
dC = 1e-3;

%%

x = 0:dx:L; x = x';
N = length(x);

%%

Sw = zeros(size(x)) + Swr; % ic
%Sw(1) = 1 - Sor; % left bc
%Sw(2) = (Sw(1) + Sw(3)) / 2;
Sw(1) = 0.64;

[~, dfw, ~, ~, lw, lo] = BrooksCorey(Sw, Swr, 1-Sor, muw, muo, lam);
lamT = kx * (lw + lo);

dpdx = -QT / A / (kx * lw(1));
p = x * dpdx - x(end) * dpdx;
vT = -lamT * dpdx;
v = vT / phi .* dfw;
dx0 = dt * v;

%%

[B, rb, rbx] = systemMatrix(x, shp, ch, mp);
Hx = mqDerivatives(rb, rbx, shp, 1);
%Hx = [Hx, ones(N, 1), zeros(N, 1)];
Dx = Hx / B(1:N, 1:N);

%%
for t = dt : dt : T
    disp(['t = ', num2str(t)]);
    
    while(true) %for ii = 1:20
        [xi, ~] = getUpstreamPoints(v, x, dx0, B, rb, rbx, shp, ch, mp, eta, dC, dt, err);
        Sw1 = interpSw(Sw, x, xi, B, rb, rbx, shp, ch, mp, eta, dC);
        
        [~, dfw, ~, ~, lw, lo] = BrooksCorey((Sw + Sw1) / 2, Swr, 1-Sor, muw, muo, lam);
        lamT = kx * (lw + lo);
        H1 = Hx;
        H1(1,:) = B(1,:);
        f = zeros(N, 1);
        f(1) = -QT / lw(1) / A / kx * lamT(1);
        vT = -B * (H1 \ f);
        v1 = vT / phi .* dfw;
        max(abs(v1-v))
        if max(abs(v1-v)) < err
            break;
        end
        v = v1;
    end
    
    Sw = Sw1;
    v = v1;
    
end
%%

