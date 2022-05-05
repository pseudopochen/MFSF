clear all; close all;
%addpath ../2D

Nx = 599;
Ny = Nx;
xv = linspace(-0.5,0.5,Nx+1);
yv = linspace(-0.5,0.5,Ny+1);
[xplt, yplt] = meshgrid(xv, yv);
uplt = zeros(size(xplt));

%Nx1 = 5199;
%Ny1 = Nx1;
%xv1 = linspace(-0.5, 0.5, Nx1 + 1);
%yv1 = linspace(-0.5, 0.5, Ny1 + 1);
%[xplt1, yplt1] = meshgrid(xv1, yv1);

t_all = [0.02, 0.06, 0.07, 0.08];

%figure('Position',[10,10,700,1000]);

[ha, pos] = tight_subplot(4,2,[.01 .03],[.01 .01],[.01 .01]);

for ii = 1:length(t_all)
  ii
  %subplot(4,2,ii*2-1); hold on;
  axes(ha(2*ii-1)); hold on;
  
  load(['dump_',num2str(t_all(ii)),'.bin']);
  uplt = interpPlotU(xplt, yplt, uplt, pnts_all, u_all, y0);
  pcolor(xplt, yplt, uplt);
  shading flat;
  axis equal; axis tight; 
  box on;
  colormap gray;
  NP = length(y0)/4;
  plot(y0(1:NP), y0(NP+1:2*NP), 'w.-','MarkerSize', 5)

  %subplot(4,2,ii*2); hold on;
  axes(ha(2*ii)); hold on;
  if ii != length(t_all)
    load(['u2_',num2str(t_all(ii)+0.01),'.bin']);
  else
    load(['u2_',num2str(0.086),'.bin']);
  endif
  pcolor(xplt, yplt, u); 
  shading flat;
  axis equal; axis tight;
  box on;
  colormap gray;
  plot(y0(1:NP), y0(NP+1:2*NP), 'w.-', 'MarkerSize', 5)  
  
endfor

%%








