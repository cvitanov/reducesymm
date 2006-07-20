function ae = kseqna(nu, a0)
%  function ae = kscvna(nu, a0) - Newton-Armijo for computing equilibria of ksfmfj

  ae = a0;  cont = 1;
  while cont,
    [g, dg] = ksfmfj(0, ae, nu);
    ng = norm(g);
%    disp(sprintf('%11.3e', ng));

    invdg = inv(dg);
    
    flag = 0;
    for ii = 0:10,      % Armijo iteration
      da = 2.^(-ii).*(invdg*g);
      a1 = ae - da;  g1 = ksfmfj(0, a1, nu);
      disp(sprintf('    %2d  %11.3e  %11.3e',ii,norm(da),norm(g1)));
      if norm(g1) < ng.*(1-2.^(-ii-1)),  ae = a1;  flag = 1;  break;  end
    end,
%    if flag == 0, break; end,
    cont = input('Continue? ');
  end,
