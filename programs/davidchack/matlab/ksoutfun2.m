function stop = ksoutfun2(a, optimValues, state)
% 2 - plot each component g_i^2 of g'*g = sum_i g_i^2

global TT A F FT NFEVAL PHASE DD HH MSTP NP AM PFLG
persistent hax hp1 hp2 hp3 hp4 hp5 hp6 hp7 nitr x k Ap ddo an han

  stop = 0;  np = 5;  danmax = 1e-2;
  switch state
  case 'init'
    k = (1:size(A,1)/2)'; an = a; han = 1;%     disp([size(k); size(A)]);
%    disp(size(A(1:2:end,end)+1i*A(2:2:end,end)));
%    disp(size(exp((-2i*pi*PHASE./DD).*k)));
    v = (PFLG.*A(1:2:end,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
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
           optimValues.funccount, TT(end), PHASE, dd(end).^2));  ddo = dd(end);
    nitr = input('NITR = ');
  case 'iter',
    v = (PFLG.*A(1:2:end-1,end)+1i*A(2:2:end,end)).*exp((-2i*pi*PHASE./DD).*k);
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
    if size(an,2)>=np, an = [an(:,2:end) a]; han = [han(2:end) optimValues.residual]; 
    else, an = [an a]; han = [han optimValues.residual]; end
    if size(an,2)>=np, dan = diff(an')';
      dand = norm(dan(:,np-1)-dan(:,np-2))./(norm(dan(:,np-1))+norm(dan(:,np-2)));
      danp = norm(dan(:,np-1)-dan(:,np-3))./(norm(dan(:,np-1))+norm(dan(:,np-3)));
    else, dand = 1.0; danp = 1.0; end
    disp(sprintf('%5d %5d %10.6f %9.6f %13.6e %13.6e %12.3e %9.1e %9.1e %9.1e', optimValues.iteration, ...
         optimValues.funccount, TT(end), PHASE, optimValues.residual, dd(end), dstp, han(end)./(han(end-1)-han(end)), dand, danp));  ddo = dd(end);
 %% curve search
    if 0 & dand < danmax & han(end)./(han(end-1)-han(end)) > 100, 
      sm = zeros(1,5); hm = sm; ni = sm;  am = zeros(length(a),5);
      [sm(1), am(:,1), hm(1), ni(1)] = crvsrch(@(am)ksfmms(am,DD,HH,MSTP,NP), an, 2);
      [sm(2), am(:,2), hm(2), ni(2)] = crvsrch(@(am)ksfmms(am,DD,HH,MSTP,NP), an, 3);
      [sm(3), am(:,3), hm(3), ni(3)] = crvsrch(@(am)ksfmms(am,DD,HH,MSTP,NP), an, 4);
      [sm(4), am(:,4), hm(4), ni(4)] = crvsrch(@(am)ksfmms(am,DD,HH,MSTP,NP), an, 5);
      [sm(5), am(:,5), hm(5), ni(5)] = crvsrch(@(am)ksfmms(am,DD,HH,MSTP,NP), an, 6);
      disp(sprintf('%12.3e %12.3e %5.1f %3d %12.3e %5.1f %3d %12.3e %5.1f %3d %12.3e %5.1f %3d %12.3e %5.1f %3d', ...
           optimValues.residual, [hm;sm;ni]));
      disp(sprintf('%12.3e %12.3e %5.1f %12.3e %5.1f %12.3e %5.1f %12.3e %5.1f %12.3e %5.1f', ...
           optimValues.residual, ...
           hm(1),(hm(1)-han(end))./(han(end)-han(end-1)), ...
           hm(2),(hm(2)-han(end))./(han(end)-han(end-1)), ...
           hm(3),(hm(3)-han(end))./(han(end)-han(end-1)), ...
           hm(4),(hm(4)-han(end))./(han(end)-han(end-1)), ...
           hm(5),(hm(5)-han(end))./(han(end)-han(end-1))));
      [hmm,ihm] = min(hm);  
      if 1 & (hm(ihm)-han(end))./(han(end)-han(end-1)) > 6,
        AM = am(:,ihm); stop = 1; end
    end
    if optimValues.iteration >= nitr,  nitr = nitr + input('NITR = ');  else,  drawnow;  end
    if optimValues.iteration == nitr,  AM = a; stop = 1; end
  case 'done', end
  