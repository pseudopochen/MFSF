[t1i, p1i] = interpShockT(t1, y1, ndec);
  
if isFirstStep
  pnts_all = [pnts_all; p1i(NPS+1:2*NPS,:)];
  [dm0, shp0] = derSpatial(pnts_all, ns, shape2, dshp, minK, maxK);
    
  t1i = t1i(2:end);
  p1i = p1i(NPS+1:end,:);
  u_all = [u_all; ul * ones(NPS,1)];
  
  isFirstStep = false;
endif

if resample
  p1i(end-NPS+1:end,1) = y0(1:ndec:NP);
  p1i(end-NPS+1:end,2) = y0(NP+1:ndec:2*NP);
endif
  
[u_all, pnts_all, dm0, shp0] = calcSat(t1i, p1i, dm0, ul, u_all, pnts_all, ns, shp0, dshp, minK, maxK, mu);

figure(1)
plot(p1i(:,1), p1i(:,2), 'k+')

figure(2); hold on; 
idx = 10;
dist = sqrt(pnts_all(idx:NPS:end,1).^2 + pnts_all(idx:NPS:end,2).^2);
plot(dist, u_all(idx:NPS:end), 'k+-')
