function status = ksdupoplot(t, a, flag);
  
global TT A F FT NFEVAL
persistent hp1 hp2 hp3 hp4 hp5 nstp nitr x

  status = 0;
  if isempty(flag),
    figure(3);
    set(hp1,'xdata',A(1,:),'ydata',A(2,:)); 
    plot(A(1,1),A(2,1),'ro',A(1,end),A(2,end),'ko');
    figure(4);
    dd = sqrt(sum((A-repmat(A(:,1),1,size(A,2))).^2));
    dn = sum(repmat(FT./norm(FT),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
    d0 = sum(repmat(F./norm(F),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
    set(hp2,'xdata',TT,'ydata', dd);  set(hp3,'xdata',TT,'ydata',dn);  set(hp4,'xdata',TT,'ydata',d0);
    plot(TT(end), dd(end), 'ko', TT(end), dn(end), 'ro', TT(end), d0(end), 'bo');
    figure(5);
    v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
    set(hp5,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    nstp = nstp + 1;
    disp(sprintf('%5d %5d %10.6f %10.6f %15.6e', nstp, NFEVAL, t, TT(end), dd(end)));
    if nstp >= nitr,  nitr = nitr + input('NITR = ');  else,  drawnow;  end
  else
    switch(flag),
    case 'init'
      NFEVAL = 0;  nstp = 1;
      figure(3); 
      clf; hp1 = plot(A(1,:),A(2,:),'.-'); grid on; hold on;
      plot(A(1,1),A(2,1),'ro',A(1,end),A(2,end),'ko');
      dd = sqrt(sum((A-repmat(A(:,1),1,size(A,2))).^2));
      dn = sum(repmat(FT./norm(FT),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
      d0 = sum(repmat(F./norm(F),1,size(A,2)).*(A-repmat(A(:,1),1,size(A,2))));
      figure(4); 
      clf; hp2 = plot(TT, dd,'.-','color',[0 .8 0]); grid on; hold on;
      hp3 = plot(TT, dn, '.-');  hp4 = plot(TT, d0, 'k.-');
      plot(TT(end), dd(end), 'ko', TT(end), dn(end), 'ro', TT(end), d0(end), 'bo');
      figure(5); 
      clf; d = 22; N = size(A,1)+2;  x = d.*(1:N)'./N;  v = A(1:2:end,:) + 1i*A(2:2:end,:);
      v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  u = real(fft(v));
      hp5 = pcolor(x,TT,u'); shading flat; axis tight;
      disp(sprintf('%5d %5d %10.6f %10.6f %15.6e', nstp, NFEVAL, t(1), TT(end), dd(end)));
      nitr = input('NITR = ');
    case 'done'
      disp('done');
    end,
  end,
