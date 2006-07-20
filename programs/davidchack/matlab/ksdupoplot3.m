function status = ksdupoplot3(t, a, flag);
% 3 - for use with ksfmfm  

global TT A F FT NFEVAL DD
persistent hp1 hp2 hp3 hp4 hp5 nstp nitr x ddo

  status = 0;
  if isempty(flag),
    F = ksfm(0, A(:,1), DD);  FT = ksfm(0, A(:,end), DD);
    figure(3);
    set(hp1,'xdata',A(1,:),'ydata',A(2,:)); 
    plot(A(1,1),A(2,1),'ro',A(1,end),A(2,end),'ko');
    figure(4);
    dd = sqrt(sum((A-repmat(A(:,1),1,size(A,2))).^2));
    dn = sum(repmat(FT./norm(FT),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
    d0 = sum(repmat(F./norm(F),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
    set(hp2,'xdata',TT,'ydata', dd);  set(hp3,'xdata',TT,'ydata',dn);  set(hp4,'xdata',TT,'ydata',d0);
    plot(TT(end),dd(end),'ko', TT(end),dn(end),'ro', TT(end),d0(end),'bo');
    figure(5);
    v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    set(hp5,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    nstp = nstp + 1;
    if ddo-dd(end) == 0, dstp = 0; else dstp = dd(end)./(ddo-dd(end)); end
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e %12.3e', nstp, ...
         NFEVAL, t, TT(end), dd(end), dstp));  ddo = dd(end);
%    if dd(end) > 1.2, status = 1; end
    if nstp >= nitr, nitr = nitr + input('NITR = '); else, drawnow; end
%    if nstp >= nitr | dstp < 0, status = 1; end
    if nstp >= nitr, status = 1; end
  else
    switch(flag),
    case 'init'
      NFEVAL = 0;  nstp = 1;
      F = ksfm(0, A(:,1), DD);  FT = ksfm(0, A(:,end), DD);
      figure(3); set(gcf,'pos',[265 525 495 420]);
      clf; hp1 = plot(A(1,:),A(2,:),'.-'); grid on; hold on;
      plot(A(1,1),A(2,1),'ro',A(1,end),A(2,end),'ko');
      dd = sqrt(sum((A-repmat(A(:,1),1,size(A,2))).^2));
      dn = sum(repmat(FT./norm(FT),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
      d0 = sum(repmat(F./norm(F),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
      figure(4); set(gcf,'pos',[760 525 635 420]);
      clf; hp2 = plot(TT, dd,'.-','color',[0 .8 0]); grid on; hold on;
      hp3 = plot(TT, dn, '.-');  hp4 = plot(TT, d0, 'k.-');
      plot(TT(end),dd(end),'ko', TT(end),dn(end),'ro', TT(end),d0(end),'bo');
      figure(5);  set(gcf,'pos',[5  35  250  420]);
      clf; N = size(A,1)+2;  x = DD.*(1:N)'./N;  v = A(1:2:end,:) + 1i*A(2:2:end,:);
      v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
      hp5 = pcolor(x,TT,u'); shading flat; axis tight;
      disp(sprintf('%5d %5d %10.6f %10.6f %15.6e', nstp, ...
           NFEVAL, t(1), TT(end), dd(end)));  ddo = dd(end);
      nitr = input('NITR = ');
    case 'done'
      disp('done');
    end,
  end,
