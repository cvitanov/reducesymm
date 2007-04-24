%% Regular and turbulent (chaotic) dynamics of Kuramoto-Sivashinsky equation

%% 
% $$ u_t = -uu_x - u_{xx} - u_{xxxx},  \quad x \in [-L/2, L/2]$$
%
% with periodic boundary condition: 
%
% $$ u(x+L,t) = u(x,t). $$
% 
% The solutions depend on the parameter _L_.
%
% Matlab function |ksfmstp.m| can be used to find numerical solution of
% the Kuramoto-Sivashinsky equation using Exponential Time Difference
% Runge-Kutta 4th order method (ETDRK4):
  N = 64;  L = 22;  h = 0.25;  nstp = 1000;
  a0 = zeros(N-2,1);  a0(1:4) = 0.6; % just some initial condition
  [tt, at] = ksfmstp(a0, L, h, nstp, 1);
  fig1 = figure('pos',[5 550 600 400],'color','w'); plot(tt,at,'.-');  
  title('Solution of Kuramoto-Sivashinsky equation with L = 22: Fourier modes');
  xlabel('Time'); ylabel('Fourier modes');
%%
% The same solution can be viewed in real space by transforming the
% solution using function |ksfm2real.m|:
  [x, ut] = ksfm2real(at, L);
  fig2 = figure('pos',[5 270 600 200],'color','w'); pcolor(tt,x,ut); shading interp;
  title('Solution of Kuramoto-Sivashinsky equation with L = 22: Real space');
  xlabel('Time'); ylabel('x','rotat',0);
  
%%
% The solution can be also viewed in the 'state subspace', spanned 
% by selected eigenvectors of one of the equilibrium points
  load('ks22statespace','v');  % v - subspace orthonormal vectors
  av = v'*at;  % projection of solution onto the subspace
  fig3 = figure('pos',[610 450 500 500],'color','w');
  plot3(av(1,:),av(2,:),av(3,:),'b.-'); grid on;
  xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0); view(-75,15);
  
%% Sensitive dependence on initial condition
% Start another solution with the nearby initial condition 
% and see how it diverges from the first one.
  a1 = a0;  a1(1) = a1(1) + 1e-3;
  [tt, at1] = ksfmstp(a1, L, h, nstp, 1);  av1 = v'*at1;
  figure(fig3); hold on;
  plot3(av1(1,:),av1(2,:),av1(3,:),'r.-');
  
%% Dependence of KS solutions on L
  LL = [10 12 16 22 24];
  fig4 = figure('pos',[500 100 600 800],'color','w');
  for ii = 1:length(LL);
    [tt, at] = ksfmstp(a0, LL(ii), h, nstp, 1);
    [x, ut] = ksfm2real(at, L);
    subplot(length(LL),1,ii);
    pcolor(tt,x,ut); shading interp;
    title(['L = ' num2str(LL(ii))]);
  end
    