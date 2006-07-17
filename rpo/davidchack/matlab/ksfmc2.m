function g = ksfmc2(t,a);
%  function g = ksfmc2(t,a);

  global NSTPS
  persistent nstpold nfeval
      
  N = 32;  nu = 0.015;  scale = 0.005;
  f = ksfm(nu, a(1:32));  %f = f./norm(f);
  opt.RelTol = 1e-10;  opt.AbsTol = 1e-10;
  opt.Jacobian = 'ksfmj_mex';  opt.JacobianLowerBandwidth = 0;  
  opt.JacobianUpperBandwidth = 0;  opt.InitialStep = 0.005;
  [s,b] = radau5Mex('ksfm_mex', [0:0.02:1]*a(33), a(1:32), opt);
%  [b,s] = flowdp8('ksfm', [32; 0.015], a(1:32), 0, [0:0.02:1]*a(33), 2e-4);
  if isempty(nfeval), nstpold = -1;  nfeval = 0; end,
  nfeval = nfeval + 1;
  ft = ksfm(nu, b(end,1:32)');%  disp(norm(ft));  ft = ft./norm(ft);
  
  g = b(end,1:32)'-a(1:32);  ng = norm(g);%  g = g./ng;
%  fg = f'*g;  
%  g = [g - fg.*f; -ft'*g];
  g = [g; -scale.*(f'*g)];
  
  fg = (f'*g(1:N))./(norm(f)*ng);
  
  %% disp([NSTPS nstpold nfeval t norm(g)]); pause;

  if NSTPS > nstpold,
    global h1 h2 h3;
    dn = sum(repmat(f./norm(f),1,size(b,1)).*(b(:,1:32)'-repmat(a(1:32),1,size(b,1))));
    dd = sqrt(sum((repmat(a(1:32),1,size(b,1))-b(:,1:32)').^2));
    figure(1);  set(h1,'xdata',s,'ydata',dn);  plot(s(end),dn(end),'ko');
    set(h2,'xdata',s,'ydata',dd); plot(s(end),dd(end),'ro');
    figure(2);  set(h3,'xdata',b(:,1),'ydata',b(:,2));
    plot(b(1,1),b(1,2),'ro',b(end,1),b(end,2),'ko');
    disp(sprintf('%5d %5d %10.6f %10.6f %10.6f', NSTPS, nfeval, t, ng, fg));
    nstpold = NSTPS;  if mod(nstpold,100) == 99; pause; else; drawnow; end,
  end,
  
  return;