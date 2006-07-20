function status = kszolplot(t, u, flag);

global NFEVAL U0 UT TT X
persistent nstp nitr hp1 hp2 hp3

  status = 0;
  if isempty(flag),
    nstp = nstp + 1; figure(2);
    set(hp1(1),'ydata',U0);  set(hp1(2),'ydata',UT);
    figure(3); set(hp2,'ydata',[get(hp2,'ydata') t],...
      'zdata',[get(hp2,'zdata'); UT'-U0'],...
      'cdata',[get(hp2,'cdata'); UT'-U0']); axis tight;
    figure(4); set(hp3,'xdata',[get(hp3,'xdata') t],'ydata',...
      [get(hp3,'ydata') norm(UT-U0)]);
    disp(sprintf('%5d %5d %10.4f %10.4f %15.4e', nstp, NFEVAL,...
      t, TT, norm(UT-U0)));
    if nstp >= nitr, nitr = nitr + input('NITR = ');  
    else, drawnow; end,
  else
    switch(flag),
    case 'init'
      NFEVAL = 0;  nstp = 1;
      figure(2); clf; hp1 = plot(X, U0, 'r.-', X, UT, 'k.-');
      set(gca,'ylim',[-3 3]); grid on;
      figure(3); clf; 
      hp2 = pcolor(X,[0 1e-6],zeros(2,length(X))); shading flat;
      figure(4);  hp3 = plot(0,norm(UT-U0),'k.-'); hold on; grid on;
      disp(sprintf('%5d %5d %10.4f %10.4f %15.4e', nstp, NFEVAL,...
        t(2), TT, norm(UT-U0)));
      nitr = 50; % nitr = input('NITR = ');
    case 'done'
      disp('done');
    end,
  end,
  