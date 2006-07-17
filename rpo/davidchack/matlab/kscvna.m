function [tp, ap] = kscvna(nu, t0, a0)
%  function [tp, ap] = kscvna(nu, t0, a0) - Newton-Armijo (Zoldi's approach)

  N = length(a0);  h = 0.002;  tp = t0;  ap = a0;  cont = 1;
  while cont,
    [tt, aa, da] = kscvfmj(nu, tp, ap, h, 0);
    g = ap - aa;  Dg = da - eye(N);  ng = norm(g);
    disp(sprintf('%8.5f  %11.3e', tp, ng));

%    f0 = ksfft(nu, ap);  ft = ksfft(nu, aa);
    [tt, aa, f0, ft] = kscvfm(nu, tp, ap, h, 0);
    invDg = inv([Dg ft; f0' 0]);
    
    flag = 0;
    for ii = 0:10,      % Armijo iteration
      da = 2.^(-ii).*(invDg*[g; 0]);
      a1 = ap + da(1:N);  te1 = tp + da(N+1);
      if te1 < 2.0*tp & te1 > 0.5*tp,
        [tt, a] = kscvfm(nu, te1, a1, h, 0);
        g1 = a1 - a;  disp(sprintf('    %2d  %11.3e  %11.3e',ii,norm(da),norm(g1)));
        if norm(g1) < ng.*(1-2.^(-ii-2)),  ap = a1;  tp = te1;  flag = 1;  break;  end
      end,
    end, 
    if flag == 0, break; end,  %cont = input('Continue? ');
  end,
