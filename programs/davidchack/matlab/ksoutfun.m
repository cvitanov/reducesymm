function stop = ksoutfun(a, optimValues, state)

global TT A F FT NFEVAL PHASE DD
persistent hp1 hp2 hp5 nitr x k Ap ddo

  stop = 0;
  switch state
  case 'init'
    k = (1:size(A,1)/2)';  
    v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
    Ap = zeros(size(A(:,1)));  Ap(1:2:end-1) = real(v);  Ap(2:2:end) = imag(v);
    figure(3); set(gcf,'pos',[265 525 495 420]);
    clf; hp1 = plot([A(1,:) Ap(1)],[A(2,:) Ap(2)],'.-'); grid on; hold on;
    plot(A(1,1),A(2,1),'ro',Ap(1),Ap(2),'ko');
    dd = sqrt(sum(([A Ap]-repmat(A(:,1),1,size(A,2)+1)).^2));
%    optimValues, disp(dd(end)^2);
    figure(4); set(gcf,'pos',[760 525 635 420]);
    clf; hp2 = plot([TT TT(end)], dd,'.-','color',[0 .8 0]); grid on; hold on;
    plot(TT(end), dd(end), 'ko');
    figure(5);  set(gcf,'pos',[5  35  250  420]);
    clf; N = size(A,1)+2;  x = DD.*(1:N)'./N;  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    hp5 = pcolor(x,TT,u'); shading flat; axis tight;
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e', optimValues.iteration, ...
           optimValues.funccount, TT(end), PHASE, dd(end)));  ddo = dd(end);
%    stop = 1;
    nitr = input('NITR = ');
  case 'iter'
    v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
    Ap(1:2:end-1) = real(v);  Ap(2:2:end) = imag(v);
    figure(3);
    set(hp1,'xdata',[A(1,:) Ap(1)],'ydata',[A(2,:) Ap(2)]); 
    plot(A(1,1),A(2,1),'ro',Ap(1),Ap(2),'ko');
    figure(4);
    dd = sqrt(sum(([A Ap]-repmat(A(:,1),1,size(A,2)+1)).^2));
    set(hp2,'xdata',[TT TT(end)],'ydata', dd);
    plot(TT(end), dd(end), 'ko');
    figure(5);
    v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    set(hp5,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    if ddo-dd(end) == 0, dstp = 0; else dstp = dd(end)./(ddo-dd(end)); end
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e %12.3e', optimValues.iteration, ...
           optimValues.funccount, TT(end), PHASE, dd(end), dstp));  ddo = dd(end);
    if optimValues.iteration >= nitr,  nitr = nitr + input('NITR = ');  else,  drawnow;  end
    if optimValues.iteration == nitr, stop = 1; end
  case 'done'
  end
  