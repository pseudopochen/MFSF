%clear all;

T = 1000; % days
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

SBL = 0.64; %0.538015;
sw0 = SBL : 0.01 : 1 - Sor; sw0 = sw0';
%[krw, kro, fw, dfw]=mobility(Sw0, Swr, Sor);
[fw, dfw, krw, kro, lm, lo] = BrooksCorey(sw0, Swr, 1-Sor, muw, muo, lam);
xSw_all = [];

figure; hold on;

for t = 600:600%dT:dT:T
    xSw = QT * t / A / phi * abs(dfw);
    xx = [max(xSw); L];
    sw1 = [Swr; Swr];
    
    xSw_all = [xSw_all, xSw];
    plot([flipud(xSw); xx] / L, [flipud(sw0); sw1], 'k-','LineWidth',1);
end

plt_ana = [[flipud(xSw); xx] / L, [flipud(sw0); sw1]];
save('-binary','dump_ana.bin','plt_ana');

%plot([0, 1], [SBL, SBL], 'k--')


