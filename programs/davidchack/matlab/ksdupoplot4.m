function status = ksdupoplot4(t, a, flag);
% 2 - with shift  
% 4 - plot each component g_i^2 of g'*g = sum_i g_i^2

global TT A F FT NFEVAL PHASE DD
persistent hax hp1 hp2 hp3 hp4 hp5 hp6 hp7 nstp nitr x k Ap ddo mdsqo flg

  status = 0;
  switch flag
  case 'init'
    NFEVAL = 0;  nstp = 1;  k = (1:size(A,1)/2)';  
    v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
    Ap = zeros(size(A(:,1))); Ap(1:2:end-1) = real(v); Ap(2:2:end) = imag(v);
    figure(3); clf; set(gcf,'pos',[5 537 1392 410]);
    hax = subplots(1,3,[0.02 0.01 0.0 0.0],[0.02 0.02 0.10 0.06]);
    axes(hax(1));  hp1 = plot([A(1,:) Ap(1)],[A(2,:) Ap(2)],'.-'); grid on; hold on;
    hp2 = plot(A(1,1),A(2,1),'ro',Ap(1),Ap(2),'ko');
    dd = sqrt(sum(([A Ap]-repmat(A(:,1),1,size(A,2)+1)).^2));
    axes(hax(2)); hp3 = plot([TT TT(end)], dd,'.-','color',[0 .8 0]); grid on; hold on;
    hp4 = plot(TT(end), dd(end), 'ko');
    axes(hax(3)); dsqi = (Ap - A(:,1)).^2; ldsq = length(dsqi); [mdsq,idsq] = max(dsqi);
    hp6 = plot(zeros(2,ldsq),[dsqi dsqi]','.-'); grid on;
    hp7 = text(zeros(ldsq,1),dsqi,num2str((1:ldsq)'),'horiz','left');
    set(hax(3),'ylim',[0 1.2*max(dsqi)]);
    figure(4); clf; set(gcf,'pos',[5  35  250  420]);
    N = size(A,1)+2;  x = DD.*(1:N)'./N;  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    hp5 = pcolor(x,TT,u'); shading flat; axis tight;
    disp(sprintf('%4d %5d %9.2e %10.6f %10.6f %12.4e %3d %12.4e', nstp, ...
         NFEVAL, t(1), TT(end), PHASE, dd(end), idsq, mdsq)); 
    ddo = dd(end);  mdsqo = mdsq;  nitr = input('NITR = ');
  case ''
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
    dsqi = (Ap - A(:,1)).^2;  [mdsq,idsq] = max(dsqi);
    for ihp = 1:length(hp6),
      hp6x = get(hp6(ihp),'xdata'); lhp6x = length(hp6x)-1;
      set(hp6(ihp),'xdata',[hp6x lhp6x],'ydata',[get(hp6(ihp),'ydata') dsqi(ihp)]);
      set(hp7(ihp),'pos',[lhp6x dsqi(ihp)]); end
    set(hax(3),'ylim',[0 1.2*max(dsqi)], 'xlim',[min(max(0,lhp6x-150),0.67*lhp6x) lhp6x]);
    figure(4);  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    set(hp5,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    nstp = nstp + 1;
    if (mdsqo-mdsq) == 0, mdstp = 0; else mdstp = mdsq./(mdsqo-mdsq); end
    if ddo-dd(end) == 0, dstp = 0; else dstp = dd(end)./(ddo-dd(end)); end
    disp(sprintf('%4d %5d %9.2e %10.6f %10.6f %10.4e %11.3e %3d %10.4e %11.3e', nstp, ...
         NFEVAL, t, TT(end), PHASE, dd(end), dstp, idsq, mdsq, mdstp));
    ddo = dd(end);  mdsqo = mdsq;
    if dd(end) > 1.2, status = 1; end
    if nstp >= nitr | (flg == 1 & dstp < 0),
      moreitr = input('NITR = ');
      if moreitr < 0, nitr = nstp - moreitr;  flg = 1;
      else, nitr = nstp + moreitr; flg = 0; end
    else, drawnow; end
    if nstp >= nitr, status = 1; end
  case 'done', end