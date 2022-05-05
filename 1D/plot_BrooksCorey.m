Swr = 0.16;
Sor = 0.2;

mu = 1e-2;
muw = mu;
muo = mu;
lam = 2;

sw = Swr : .01 : 1-Sor;

[fw, dfw, krw, kro, lm, lo] = BrooksCorey(sw, Swr, 1-Sor, muw, muo, lam);

figure;

[ax,h1,h2] = plotyy(sw,fw,sw,dfw);

set(ax,'YColor','k');
set(h1,'color','k','LineWidth',1);
set(h2,'color','k','LineStyle','--','LineWidth',1);
xlabel('water saturation','FontSize',18)
grid on; box on;
%print -dpng BrooksCorey.png
print -deps BrooksCorey.eps
