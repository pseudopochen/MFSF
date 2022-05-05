clear all;

T = 100; % days
dT = 100; % save result per dT days
dx = 0.125; % m
dt = 0.5; % day
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

dpdx = -QT / A / kx / lw(1);
p = x * dpdx - x(end) * dpdx;

a = kx * dt / dx^2 ;

%%

Sw0 = Sw;
Sw1 = Sw;
p0 = p;
p1 = p;

Sw_all = [];
isave = 0;
for t = dt : dt : T
  
  disp(['t = ', num2str(t)]);

  icount = 0;
  do
    Sw = Sw1;
    p = p1;
  
    dSw = a / phi * (lw(2:end-1).*p1(3:end) - (lw(2:end-1) + lw(1:end-2)) .* p1(2:end-1) + lw(1:end-2).*p1(1:end-2));
    Sw1(2:end-1) = dSw + Sw0(2:end-1);
    Sw1(1) = 1 - Sor;
    Sw1(end) = Swr;

    [~, dfw, ~, ~, lw, lo] = BrooksCorey((Sw0 + Sw1) / 2.0, Swr, 1-Sor, muw, muo, lam);
    ll = lw + lo;
  
    D = zeros(length(x), length(x));
    for ii = 2:length(x) -1
      D(ii, ii - 1) = ll(ii - 1);
      D(ii, ii + 1) = ll(ii);
      D(ii, ii) = -ll(ii - 1) - ll(ii);
    end
    D = D * a;
    D(1, 1) = -1;
    D(1, 2) = 1;
    D(end, end) = 1;
  
    b = zeros(size(x));
    b(1) = -QT / lw(1) / A / kx * dx;
    %b(end) = -QT / ll(1) / A / kx * dx;
  
    p1 = D \ b;
    %pause

    icount = icount + 1;
    disp(['icount = ', num2str(icount), ' ', num2str(max(abs(Sw-Sw1))), ...
	  ' ', num2str(max(abs(p-p1)))])
  until (max(abs(Sw - Sw1)) < err && max(abs(p - p1)) < err)

  Sw0 = Sw1;
  p0 = p1;
  
  if mod(t, dT) == 0
    Sw_all = [Sw_all, Sw1];
    isave = isave + 1;
  end
  
end

%%

%for jj = 1:isave
%  plot(x / L, Sw_all(:,jj), 'k')
%end

%save BL_analytic_fd.ou Sw_all -ascii -double

