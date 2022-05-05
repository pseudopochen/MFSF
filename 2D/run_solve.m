t = 0;
FinalTime = 8.6e-2;
isFirstStep = true;
resample = false;

y0 = [pnts(1,:)'; pnts(2,:)'; nn(1,:)'; nn(2,:)'];
%shape3 = 30;
thresh = 0.1;

pnts_all = pnts(:, 1:ndec:end)';
u_all = ones(NP/ndec,1);

while t < FinalTime
  
  dt = 1e-2;
  %if t >= 0.07
  %  dt = 5e-3
  %endif
  
  if t + dt > FinalTime
    dt = FinalTime - t;
  end
  
  disp(['t = ', num2str(t), ' to ', num2str(t+dt)]);
  
  % solve shock
  
  if t < 0.06
    nbs = neighbors_shock;
  else
    nbs = neighbors_shock_c;
  endif
  [t1, y1] = ode45(@(t,y)dydt(t,y, c, nbs, shape1, nl, useCorners, corners, dshp, minK, maxK*10), [0, dt], y0);  
  
  y0 = y1(end,:)';
  %if t >= 0.06
  [y0, resample] = adaptShockPnts(y0, isClosed, useCorners, corners, thresh);
  %endif
  
  figure(1);
  plot(y0(1:NP), y0(NP+1:2*NP), 'r.-', 'MarkerSize',10);
  
  % solve saturation
  
  solve_sat; 
  
  %% save for plot
  
  fnam = ['dump_', num2str(t), '.bin'];
  save('-binary', fnam, 'y0', 'pnts_all', 'u_all');
  
  %% increase time
  
  t = t + dt;
  
endwhile
