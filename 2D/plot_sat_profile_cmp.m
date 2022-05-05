Nx = 599;
Ny = Nx;
xv = linspace(-0.5,0.5,Nx+1);
yv = linspace(-0.5,0.5,Ny+1);
[xplt, yplt] = meshgrid(xv, yv);
uplt = zeros(size(xplt));

Nx1 = 5199;
Ny1 = Nx1;
xv1 = linspace(-0.5, 0.5, Nx1 + 1);
yv1 = linspace(-0.5, 0.5, Ny1 + 1);
[xplt1, yplt1] = meshgrid(xv1, yv1);

%%

load(['dump_',num2str(0.05),'.bin']);
uplt = interpPlotU(xplt, yplt, uplt, pnts_all, u_all, y0);

load(['u2_',num2str(0.06),'.bin']);

xi = 0:.0001:.5;
yi = xi;
dd = sqrt(xi.^2 + yi.^2);

ui = interp2(xplt, yplt, uplt, xi, yi);
u2i = interp2(xplt, yplt, u, xi, yi);

figure; hold on;
plot(dd, u2i, 'k:', 'LineWidth',3)
plot(dd, ui, 'k-', 'LineWidth',2)
xlim([0 0.7]); ylim([-0.1 1.1]); box on; grid on;
xlabel('Distance', 'FontSize', 14);
ylabel('Saturation', 'FontSize',14);
lgd = legend('FV', 'MFSF');
legend('boxoff')
set(lgd, 'fontsize',16);
print -deps sat_profile_cmp.eps
%print -dpng sat_profile_cmp.png
