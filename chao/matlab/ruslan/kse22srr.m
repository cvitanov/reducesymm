%% Symmetry-Reduced Representation (SRR) of KSE solutions

%% Kuramoto-Sivashinsky Equation (KSE)
% $$ u_t = -uu_x - u_{xx} - u_{xxxx},  \quad x \in [-L/2, L/2]$$
%
% with periodic boundary condition: 
%
% $$ u(x+L,t) = u(x,t). $$
% 
% In Fourier representation
%
% $$ u(x,t)=\sum_{k=-\infty}^{+\infty} a_k (t) e^{ i q_k x } $$
%
% where
%
% $$ q_k = 2\pi k/L $$
%
% the KSE takes the form
%
% $$ \dot{a}_k = v_k(a) = (q_k^2 - q_k^4)a_k - \frac{iq}{2}\sum_{m=-\infty}^{\infty}
% a_m a_{k-m}\,,\quad a_k \in {\cal C} $$
%
% It is convenient to represent complex modes as pairs of real variables, either
% in cartesian or polar coordinates:
%
% $$ a_k = (b_k, c_k) = b_k + ic_k = (r_k, \theta_k) = r_k e^{i\theta_k}\,. $$
%
% *Symmetries:*  
%
% If 
%
% $$ u(x,t) $$
% 
% is a solution of the KSE, then so are
% 
% $$ \tau_{\ell/L} u(x,t) = u(x+\ell,t) \quad\mbox{and}\quad R\,u(x,t) = -u(-x,t). $$
%
% The action of symmetry transformations on Fourier modes is as follows:
%
% $$ \tau_{\ell/L} a_k = e^{iq_k \ell} a_k = (r_k, \theta_k + q_k \ell) \,,\qquad R\,a_k = -a_k^\ast = (-b_k, c_k)\,. $$ 

%% Equilibria and Travelling Waves
% Equilibria and Travelling Waves of KSE for L = 22 have been discussed 
% <http://www.maths.le.ac.uk/people/rld8/temp/kse22explore.html here>.

%% Fundamental domain for KSE solutions
%
% The fundamental domain (FD) of KSE solutions is defined in Fourier space as
% follows:
%
% $$\mbox{If~~} r_1 > 0 \mbox{~~then~~} \theta_1 = \pi/2 
%   \mbox{~~and~~} \dot{\theta}_1 \leq 0\,.$$
%
% $$\mbox{If~~} r_1 = 0 \mbox{~~then~~} \arg \dot{a}_1 = \pi/2\,.$$
%
% *Note 1:* The invariant subspace of antisymmetric solutions 
%
% $$ u(x,t) = R\,u(x,t) = -u(-x,t) \quad\Rightarrow\quad Re\,a_k = 0 $$
%
% belongs to the FD. 
%
% *Note 2:* If both 
%
% $$ a_1 = 0 \mbox{~~and~~} \dot{a}_1 = 0 $$
%
% then the solution lives in the L/2-periodic invariant subspace, and the FD
% is defined in terms of 
%
% $$ r_2 \mbox{~~and~~} \theta_2\,. \mbox{~~And so on...} $$
%
% A KSE solution
%
% $$ (r_1, \theta_1, r_2, \theta_2, \ldots) $$
%
% is mapped into the FD by translation 
% 
% $$ \ell/L = \frac{\pi/2 - \theta_1}{2\pi} $$
%
% and reflection if
%
% $$ \dot{\theta}_1 > 0\,, \mbox{~~or~~} b_1\dot{c}_1 - c_1\dot{b}_1 > 0\,, \mbox{~~since} $$
%
% $$ \dot{\theta_k} = \frac{b_k\dot{c}_k - c_k\dot{b}_k}{b_k^2 + c_k^2} $$
%
% The map into FD is characterized by two parameters: the translation
% parameter
% 
% $$ \ell = \frac{\pi/2 - \theta_1}{2\pi}L $$
%
% and the reflection parameter
% 
% $$ \rho = 1 \mbox{~~when~~} \dot{\theta}_1 \leq 0 \mbox{~~and~~} -1
% \mbox{~~when~~} \dot{\theta}_1 > 0 $$
%
% Here's an example of the KSE solution and its image in the fundamental
% domain.

%%
  clear; L = 22; N = 32; h = 0.25; d = L;
  a0 = zeros(N-2,1); randn('seed',22010004); a0(1:8) = 0.2*randn(8,1);
  [tt, aa] = ksfmstp2(a0, L, h, 400, 1);  ii = 82:length(tt);
  set(gcf,'pos',[350 650 700 300]); clf; 
  plot(tt(ii)-tt(ii(1)),aa(1:6,ii),'.-'); set(gca,'pos',[0.05 0.15 0.92 0.75],'xlim',[0 95]);
  legend('b_1','c_1','b_2','c_2','b_3','c_3');
  title('KSE solution: first three modes as functions of time'); xlabel('t');
%%
  va = aa(1:2:end,ii) + 1i*aa(2:2:end,ii); theta1 = angle(va(1,:));
  fa = ksfm2(aa(:,ii),L);  thetadot1 = (aa(1,ii).*fa(2,:)-aa(2,ii).*fa(1,:))./(aa(1,ii).^2+aa(2,ii).^2);
  rho = ones(1,size(fa,2)); rho(thetadot1 > 0) = -1;
  ell = mod((pi/2 - theta1)*(L/2/pi)+L/2,L)-L/2;
  clf; plot(tt(ii)-tt(ii(1)),[ell; rho]','.-'); 
  set(gca,'pos',[0.05 0.15 0.92 0.75],'xlim',[0 95],'ylim',[-L L]/2); grid on;
  lh = legend('$\ell$','$\rho$'); set(lh,'interp','latex');
  title('Variables required to map the solution into the fundamental domain'); xlabel('t');
%%
  tau = exp(1i*(1:N/2-1)'*(pi/2-theta1));  vf = tau.*va;  
  af = zeros(N-2,size(vf,2));  af(1:2:end,:) = real(vf);  af(2:2:end,:) = imag(vf);
  ir = find(rho == -1);  af(1:2:end,ir) = -af(1:2:end,ir);
  plot(tt(ii)-tt(ii(1)),af(1:6,:),'.-'); set(gca,'pos',[0.05 0.15 0.92 0.75],'xlim',[0 95]);
  legend('b_1','c_1','b_2','c_2','b_3','c_3');
  title('Image of the KSE solution in the fundamental domain'); xlabel('t');
%%
% This looks rather messy, but it allows us to map all solutions related by
% the symmetries of the KSE into a single solution in the FD.

%% Images of RPOs and PPOs in the fundamental domain
% The following examples illustrate how RPOs and PPOs appear in the FD.
%%
  clear; load ks22f90h25t100; h = 0.1;
  irpo = 5;  a0 = rpo(irpo).a1; tend = 2*rpo(irpo).T1;
  [tt,aa] = ksfmetd2(a0, L, h, tend, 1);
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(['RPO: $T = $~' num2str(tend/2) '~~$\ell = $~' num2str(rpo(irpo).s1)],'interp','latex'); grid on;
  [af, ell, rho] = ksfm2fd(aa, L);
 axes(hax(2)); plot(tt,[ell; rho]','.-'); 
  set(gca,'xlim',[0 tend],'ylim',[-L L]/2,'xtick',(0:2)*tend/2); grid on;
 axes(hax(3)); plot(tt,af(1:6,:),'.-'); 
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2); grid on;
%%
% The RPO is shown over two periods.  In the FD it looks like a periodic
% orbit.  The translation variable correctly trances the shift of the RPO, 
% and so the difference between its values at the beginning and at the end of 
% the period is equal to the RPO shift (with a minus sign).
%%
  ippo = 6;  a0 = ppo(ippo).a1; tend = 2*ppo(ippo).T1;
  [tt,aa] = ksfmetd2(a0, L, h, tend, 1);
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(sprintf('PPO: $T = $%8.4f',tend/2),'interp','latex'); grid on;
  [af, ell, rho] = ksfm2fd(aa, L);
 axes(hax(2)); plot(tt,[ell; rho]','.-'); 
  set(gca,'xlim',[0 tend],'ylim',[-L L]/2,'xtick',(0:2)*tend/2); grid on;
 axes(hax(3)); plot(tt,af(1:6,:),'.-'); 
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2); grid on;
%%
% In case of the PPO, it is clear that it also becomes a periodic orbit in
% FD and the translation variable shows the 'winding' and 'unwinding' of the
% shift over two periods.

%% Poincare map of the KSE solutions in the fundamental domain
% *Idea:* Maybe we can define the Poincare map in the FD using points
% along the KSE flow where the reflection parameter changes (both from 1
% to -1 and from -1 to 1).  Obviously, this happens whenever
%
% $$ \dot{\theta}_1 = 0\,. $$
%
% Note that RPOs and PPOs will be represented by periodic orbits of this
% Poincare map.  Besides, all RPOs will be represented by periodic orbits with
% only _even_ periods while PPOs will be represented only by _odd_ period
% orbits.

%% Using equilibria as reference states for constructing SRR
% *Idea:* Instead of using the Fourier mode to quotient the SO(2) symmetry,
% use projections onto the equilibria E1, E2, and/or E3 to define and
% subtract the shift.  
%
% In order to see how this can be done, it is convenient to work with
% complex-valued vectors
%
% $$ a = \{a_k, k = 1,\ldots,N\} \in {\cal C}^N $$
%
% $$ \theta_1 = \arg (e_1^\dagger\, a) $$
%
% where
%
% $$ e_1 = (1+0i, 0, 0, \ldots) $$
%
% is the basis vector corresponding to the first Fourier mode and 'dagger'
% denotes Hermitian transpose.
%
% In the same way we can define a shift with respect to any state of 
% the KSE, e.g. one of the equilibrium states:
%
% $$ \theta_q = \arg (a_q^\dagger\,a) $$
%

%% How to quotient SO(2) symmetry
% In order to quotient the SO(2) symmetry we need to be able to define,
% for any state of the KS system _u(x)_ (_a_ in Fourier space), 
% a 'shift' parameter
%
% $$ s(a) \in S^1 $$
%
% such that
%
% $$ \tau_{-s(a)} a $$
%
% is in M/SO(2).
%
%   This
% parameter must satisfy the following _monotonicity condition_ with 
% respect to the translation of the state _a_:
%
% $$ s(\tau_{\ell/L}a) = s(a) + \phi(\ell/L) $$
%
% where
%
% $$ \phi(x) $$
% 
% should be a continuous strictly monotonic function on the domain
%
% $$ x \in [0, 1)\,. $$
%
% This condition is necessary to avoid any ambiguity in the definition of
% the shift parameter.
%
% This condition is clearly satisfied by the first Fourier mode, for which
% 
% $$ \arg (e_1^\dagger\,\tau_{\ell/L}a) = \theta_1 + 2\pi\ell/L\,. $$
%
% It is also possible to produce a unique shift parameter by the difference
% between phases of modes k and k+1:
%
% $$ \arg (e_{k+1}^\dagger\,\tau_{\ell/L}a) - 
% \arg (e_{k}^\dagger\,\tau_{\ell/L}a) = \theta_{k+1} + 2\pi (k+1)\ell/L 
% - \theta_k - 2\pi k\ell/L = \,. $$
%
% At the same time, when a reference state mixes different FMs
% (e.g. the equilibrium states), I don't believe it is possible to
% satisfy the monotonicity condition.

%%
  load ks22uss2a;  vref1 = a0(1:2:end) + 1i*a0(2:2:end);
  load ks22uss3a;  vref2 = a0(1:2:end) + 1i*a0(2:2:end);
  L = 22; N = 32; h = 0.05; d = L;
%  vref1 = zeros(15,1) + 1i*zeros(15,1); vref1(2) = 0 + 1i;
%  vref2 = zeros(15,1) + 1i*zeros(15,1); vref2(3) = 0 + 1i;
  a0 = zeros(N-2,1); randn('seed',22010005); a0(1:8) = 0.2*randn(8,1);
  [tt, aa] = ksfmstp2(a0, L, h, 2000, 1);  ii = 410:length(tt);
  tt = tt(ii)-tt(ii(1));  aa = aa(:,ii);  [x, uu] = ksfm2real(aa, L, 64);
  set(gcf,'pos',[350 350 700 600]); clf; 
%  subplot(1,4,1); pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
  va = aa(1:2:end,:) + 1i*aa(2:2:end,:);
%  subplot(1,2,2);
  for shift = 0:0.05:1, tau = exp(2i*pi*(1:15)'*shift); 
    phase = mod(angle((tau.*vref1)'*va)-angle((tau.*vref2)'*va)+pi,2*pi)-pi;
    plot(phase,tt,'.'); hold all; end

%%  
  clear; L = 22; N = 32; h = 0.25; d = L;
  a0 = zeros(N-2,1); randn('seed',22010005); a0(1:8) = 0.2*randn(8,1);
  load ks22uqo033b;
  [tt, aa] = ksfmetd2(a0, L, h, tend, 1);  ii = 1:length(tt);
  tt = tt(ii)-tt(ii(1));  aa = aa(:,ii);  [x, uu] = ksfm2real(aa, L, 64);
  va = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  set(gcf,'pos',[250 350 800 600]); clf; 
  hax = subplots(1,4,[0 0 0 0],[0.02 0.02 0.08 0.04]);
 axes(hax(1));
  tau1 = exp(2i*pi*(.1:.1:1)'); 
  phase1 = mod(angle(tau1*va(1,:))+pi,2*pi)-pi;
  plot(phase1,tt,'.'); set(gca,'xlim',[-pi pi]);
 axes(hax(2));
  tau2 = exp(4i*pi*(.1:.1:1)'); 
  phase2 = mod(angle(tau2*va(2,:))-angle(tau1*va(1,:))+pi,2*pi)-pi;
%  phase2 = mod(angle(tau2*va(2,:))+pi,2*pi)-pi;
  plot(phase2,tt,'.'); set(gca,'xlim',[-pi pi]);
 axes(hax(3));
  tau3 = exp(6i*pi*(.1:.1:1)'); 
  phase3 = mod(angle(tau3*va(3,:))-angle(tau2*va(2,:))+pi,2*pi)-pi;
%  phase3 = mod(angle(tau3*va(3,:))+pi,2*pi)-pi;
  plot(phase3,tt,'.'); set(gca,'xlim',[-pi pi]);
 axes(hax(4));
  tau4 = exp(8i*pi*(.1:.1:1)'); 
  phase4 = mod(angle(tau4*va(4,:))-angle(tau3*va(3,:))+pi,2*pi)-pi;
%  phase4 = mod(angle(tau4*va(4,:))+pi,2*pi)-pi;
  plot(phase4,tt,'.'); set(gca,'xlim',[-pi pi]);
