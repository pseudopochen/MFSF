% following the example in 
% "An iterative adaptive multiquadric radial basis function method for the 
% detection of local jump discontinuities" 
% by Vincent R. Durante, Jae-Hun Jung (2008);

%clear all
addpath ~/Downloads/mqMatlabAll

%%

N = 50;
NR = 500;
shp = 2.6;
ch = 2;
mp = 2;

%%

x = linspace(-1, 1, N); x = x';
y = -1 * (x >= -1 & x <= -0.7) + 0 * (x > -0.7 & x <= -0.3) + 5 * (x > -0.3 & x <= 0) + 1 * (x > 0 & x <= 0.6) -1 * (x > 0.6 & x <= 0.8) + 15 * (x > 0.8 & x <= 1);

[B, rb, rbx] = systemMatrix(x, shp, ch, mp);

%%

xx = linspace(-1, 1, NR);
[yy1, alpha3, C01, C1] = interpRPIM(xx, x, y, shp, ch, mp, B, rb, rbx, 0.5, 0.001);

%%

subplot(2,1,1); cla; hold on
plot(x, y, 'k');
plot(xx, yy1, 'r')
subplot(2,1,2); cla; hold on
stem(x, C01, 'k')
stem(x, C1, 'r')