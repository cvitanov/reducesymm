function stop = ksoutfun2(a, optimValues, state)
% 2 - plot each component g_i^2 of g'*g = sum_i g_i^2

global TT A F FT NFEVAL PHASE DD
persistent hax hp1 hp2 hp3 hp4 hp5 hp6 hp7 nitr x k Ap ddo

  stop = 0;
  switch state
  case 'init'
    k = (1:size(A,1)/2)';  
    v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
    Ap = zeros(size(A(:,1)));  Ap(1:2:end-1) = real(v);  Ap(2:2:end) = imag(v);
    figure(3); clf; set(gcf,'pos',[5 537 1392 410]);
    hax = subplots(1,3,[0.02 0.01 0.0 0.0],[0.02 0.02 0.10 0.06]);
    axes(hax(1)); hp1 = plot([A(1,:) Ap(1)],[A(2,:) Ap(2)],'.-'); grid on; hold on;
    hp2 = plot(A(1,1),A(2,1),'ro',Ap(1),Ap(2),'ko');
    dd = sqrt(sum(([A Ap]-repmat(A(:,1),1,size(A,2)+1)).^2));
    axes(hax(2)); hp3 = plot([TT TT(end)], dd,'.-','color',[0 .8 0]); grid on; hold on;
    hp4 = plot(TT(end), dd(end), 'ko');
    axes(hax(3)); dsqi = (Ap - A(:,1)).^2; ldsq = length(dsqi);
    hp6 = plot(zeros(2,ldsq),[dsqi dsqi]','.-'); grid on;
    hp7 = text(zeros(ldsq,1),dsqi,num2str((1:ldsq)'),'horiz','left');
    set(hax(3),'ylim',[0 1.2*max(dsqi)]);
    figure(4); clf; set(gcf,'pos',[5  35  250  420]);
    N = size(A,1)+2;  x = DD.*(1:N)'./N;  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    hp5 = pcolor(x,TT,u'); shading flat; axis tight;
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e', optimValues.iteration, ...
           optimValues.funccount, TT(end), PHASE, dd(end)));  ddo = dd(end);
    nitr = input('NITR = ');
  case 'iter'
    v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
    Ap(1:2:end-1) = real(v);  Ap(2:2:end) = imag(v);
    set(hp1,'xdata',[A(1,:) Ap(1)],'ydata',[A(2,:) Ap(2)]);
    set(hp2(1),'xdata',[get(hp2(1),'xdata') A(1,1)],...
               'ydata',[get(hp2(1),'ydata') A(2,1)]);
    set(hp2(2),'xdata',[get(hp2(2),'xdata') Ap(1)],...
               'ydata',[get(hp2(2),'ydata') Ap(2)]);    
    dd = sqrt(sum(([A Ap]-repmat(A(:,1),1,size(A,2)+1)).^2));
    set(hp3,'xdata',[TT TT(end)],'ydata', dd);
    set(hp4,'xdata',[get(hp4,'xdata') TT(end)],'ydata',[get(hp4,'ydata') dd(end)]);
    dsqi = (Ap - A(:,1)).^2;
    for ihp = 1:length(hp6),
      hp6x = get(hp6(ihp),'xdata'); lhp6x = length(hp6x)-1;
      set(hp6(ihp),'xdata',[hp6x lhp6x],'ydata',[get(hp6(ihp),'ydata') dsqi(ihp)]);
      set(hp7(ihp),'pos',[lhp6x dsqi(ihp)]); end
    set(hax(3),'ylim',[0 1.2*max(dsqi)], 'xlim',[min(max(0,lhp6x-150),0.67*lhp6x) lhp6x]);
    figure(4);  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    set(hp5,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    if ddo-dd(end) == 0, dstp = 0; else dstp = dd(end)./(ddo-dd(end)); end
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e %12.3e', optimValues.iteration, ...
           optimValues.funccount, TT(end), PHASE, dd(end), dstp));  ddo = dd(end);
    if optimValues.iteration >= nitr,  nitr = nitr + input('NITR = ');  else,  drawnow;  end
    if optimValues.iteration == nitr, stop = 1; end
  case 'done', end
  