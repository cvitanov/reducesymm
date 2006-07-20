function g = ksfmc4(t,a);
%  function g = ksfmc4(t,a);

  global NSTPS C
  persistent nstpold nfeval nitr
  
  N = 32;  nu = 0.015;  alpha = 0.01;  h = 2e-4;
  f = ksfm(nu, a(1:N));%  f = f./norm(f);
  s = (0:0.01:1)'*a(N+1);
  b = ksfmflowmapj(nu, s, a(1:N), h);
  if isempty(nfeval), nstpold = -1;  nfeval = 0; end,
  nfeval = nfeval + 1;
  ft = ksfm(nu, b(:,end));% ft = ft./norm(ft);
  
  g = b(:,end)-a(1:N);  ng = norm(g);
  g = [C*g; -alpha.*(ft'*g)];
  
%  fg = (ft'*g(1:N))./(ng.*norm(ft));
  
  if isempty(nitr), nitr = 0; end,
  
  %% disp([NSTPS nstpold nfeval t norm(g)]); pause;
  % NSTPS = NSTPS + 1;
  if NSTPS > nstpold,
    global h1 h2 h3;
    dn = sum(repmat(ft./norm(ft),1,size(b,2)).*(b(1:N,:)-repmat(a(1:N),1,size(b,2))));
    dd = sqrt(sum((repmat(a(1:N),1,size(b,2))-b(1:N,:)).^2));
    figure(1);  set(h1,'xdata',s,'ydata',dn);  plot(s(end),dn(end),'ro');
    set(h2,'xdata',s,'ydata',dd); plot(s(end),dd(end),'ko');
    figure(2);  set(h3,'xdata',b(1,:),'ydata',b(2,:));
    plot(b(1,1),b(2,1),'ro',b(1,end),b(2,end),'ko');
    disp(sprintf('%5d %5d %10.6f %10.6f %10.6f', NSTPS, nfeval, t, ng, a(N+1)));
    nstpold = NSTPS;  
    if nstpold > nitr,  nitr = nitr + input('NITR = ');  else,  drawnow;  end,
  end,
  
  return;