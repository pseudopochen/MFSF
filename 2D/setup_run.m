addpath ~/Desktop/mqMatlabAll/
addpath ../2D/

%pkg load matgeom
warning('off');

clear all; %close all;

%% set up shock 

R = 0.02;
dtheta = 1.25;  

[pnts, tn, nn] = circleTangentNormal (R, 0:dtheta:360-dtheta);
[~,NP] = size(pnts);

figure(1); hold on;
plot(pnts(1,:),pnts(2,:), '.-'); axis equal
xlim([-0.5 0.5]); ylim([-0.5 0.5]); grid on; box on;

##load('../2D/y_all.bin','-binary');
##figure(1); hold on; 
##for ii = 1:9 
##  plot(y_all(1:72,ii), y_all(72+1:2*72,ii), 'k.-');
##end


% shock speed init

ul = 0.5773503;
mu = 0.5;
[fl, dfl] = fractionalFlow(ul, mu);
c = fl / ul;

w = shockSpeed(pnts, nn, c);

% shock point neighbors

corners = 1:(45/dtheta):NP; %[1, 10, 19, 28, 37, 46, 55, 64];
useCorners = true;
isClosed = true;

nl = 2;
nr = 2;

neighbors_shock = [];
for ii = 1:NP
  [idxs, iroi, ~] = getSupportDomain(ii, nl, nr, NP, isClosed, false, corners);
  neighbors_shock = [neighbors_shock; idxs];
end

neighbors_shock_c = [];
for ii = 1:NP
  [idxs, iroi, ~] = getSupportDomain(ii, nl, nr, NP, isClosed, true, corners(2:2:end));
  neighbors_shock_c = [neighbors_shock_c; idxs];
end

shape1 = 1e-2;

%% debug shock setup 

%setup_debug;

%% set up saturation

ns = 11;
minK = 1e+12;  
maxK = 1e+15;  
dshp = 0.05;
shape2 = 21.3;
ndec = round(5/dtheta);
NPS = NP / ndec;

%% set up plotting grid

Nx = 599;
Ny = Nx;
xv = linspace(-0.5,0.5,Nx+1);
yv = linspace(-0.5,0.5,Ny+1);
[xplt, yplt] = meshgrid(xv, yv);
uplt = zeros(size(xplt));

%% time stepping
run_solve;


