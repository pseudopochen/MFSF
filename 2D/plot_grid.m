Bnd = [-0.5, -0.5; 0.5, 0.5];
%[Xn, bpi] = nearOptimalCenters(50,50,250,0.02,5,Bnd,30,30);
load ./Xn.bin;

t_all = [0.02, 0.06, 0.07, 0.08];

[ha, pos] = tight_subplot(2,2,[.01 .03],[.01 .01],[.01 .01]);

for ii = 1:length(t_all)
  axes(ha(ii)); hold on;
  
  fnam = ['dump_',num2str(t_all(ii)),'.bin']
  load(fnam,'-binary');
  NP = length(y0)/4;
  
  umask = inpoly2(Xn, [y0(1:NP),y0(NP+1:2*NP)]);
  Xnout = Xn(!umask,:);
  
  plot(Xnout(:,1),Xnout(:,2),'k.','MarkerSize',10);
  plot(pnts_all(:,1),pnts_all(:,2),'k.','MarkerSize',10);
  
  plot(y0(1:NP), y0(NP+1:2*NP), 'k.-');
  axis equal; axis tight; box on; grid on;
endfor
