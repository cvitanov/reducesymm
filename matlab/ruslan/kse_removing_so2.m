%% Mapping KSE solutions to M/SO(2)

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


%% Mapping KSE solutions to M/SO(2) - Ruslan's way
% We can map a KSE solution to M/SO(2) using the phase of the 1st Fourier 
% mode.  That is, a KSE solution
%
% $$ a = (b_1, c_1, b_2, c_2, \ldots) = (r_1, \theta_1, r_2, \theta_2, \ldots) $$
%
% is mapped into M/SO(2) by translation 
% 
% $$ \tau_{\ell/L},\quad \mbox{where}\quad \ell/L = \frac{\pi/2 - \theta_1}{2\pi} $$
%
% So the phase of the 1st Fourier mode of the mapped solution remains constant
%
% $$ \tau_{\ell/L} a =  (0, \bar{c}_1, \bar{b}_2, \bar{c}_2, \ldots) = (r_1, \pi/2, r_2, \bar{\theta}_2, \ldots) $$
%

%% Images of RPOs and PPOs in the fundamental domain
% The following examples illustrate how RPOs and PPOs appear in the M/SO(2).
%%
  clear; load ks22f90h25t100; h = 0.25;
  irpo = 5;  a0 = rpo(irpo).a; tend = 2*rpo(irpo).T;
  [tt,aa] = ksfmetd2(a0, L, h, tend, 1);
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(['RPO: $T = $~' num2str(tend/2) '~~$\ell = $~' num2str(rpo(irpo).s)],'interp','latex'); grid on;
%  legend('b_1','c_1','b_2','c_2','b_3','c_3','location','eastoutside');
  [af, ell] = ksfm_so2(aa, L);
 axes(hax(2)); plot(tt, ell*L,'.-'); 
  set(gca,'xlim',[0 tend],'ylim',[-L L]/2,'xtick',(0:2)*tend/2); grid on;
 axes(hax(3)); plot(tt,af(1:6,:),'.-'); 
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2); grid on;
%%
% The RPO is shown over two periods.  In M/SO(2) it looks like a periodic
% orbit.  The shift variable correctly traces the shift of the RPO, 
% and so the difference between its values at the beginning and at the end of 
% the period is equal to the RPO shift (with a minus sign).
%%
  ippo = 6;  a0 = ppo(ippo).a; tend = 2*ppo(ippo).T;
  [tt,aa] = ksfmetd2(a0, L, h, tend, 1);
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(sprintf('PPO: $T = $%8.4f',tend/2),'interp','latex'); grid on;
  [af, ell] = ksfm_so2(aa, L);
 axes(hax(2)); plot(tt,ell*L,'.-'); 
  set(gca,'xlim',[0 tend],'ylim',[-L L]/2,'xtick',(0:2)*tend/2); grid on;
 axes(hax(3)); plot(tt,af(1:6,:),'.-'); 
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2); grid on;
%%
% In the case of a PPO, it is clear that it also becomes a periodic orbit in
% M/SO(2) and the shift variable shows the 'winding' and 'unwinding' of the
% shift over two periods.

%% Mapping KSE solutions to M/SO(2) - Predrag's way
% We can map the solution of KSE starting from a_0 to M/SO(2) by
% removing from the flow _v(a)_ component tangent to the SO(2) group 
% action at a_0
%
% $$ t(a) = (t_1, t_2, \ldots)^\mathsf{T}, \quad t_k = ik(a_0)_k $$
%
% In the representation
% 
% $$ a_0 = (b_1, c_1, b_2, c_2, \ldots)^\mathsf{T},\quad t(a) = G a_0, \quad\mbox{where} $$
%
% $$ G = \left(\begin{array}{ccccc} 0 & 1 & 0 & 0 & \cdots\\ -1 & 0 & 0 & 0 & \cdots\\0 & 0 & 0 & 2 & \cdots\\0 & 0 & -2 & 0 & \cdots \\ \vdots & \vdots & \vdots & \vdots & \ddots \end{array}\right) $$
%
% The component of _v(a)_ tangent to _t(a_0)_ is given by
%
% $$ v_{\|}(a) = P_{\|}(a_0) v(a)\,\quad\mbox{where} P_{\|}(a) = \frac{t\,t^\mathsf{T}}{t^\mathsf{T} t} $$
%
% The KSE solution staying within M/SO(2) can be obtained by solving
% the modified equation
%
% $$ \dot{a}_\perp = v(a) -  v_{\|}(a)$$
%
% Here's an example of this approach applied to the same RPO and PPO.  
% To keep the calculations as explicit as possible, I'm going to use
% Euler's method to solve the original and the modified equations.
%
%%  
  clear; load ks22f90h25t100;  L = 22;  N = 32;  irpo = 5;  h = 0.05;
  a0 = rpo(irpo).a;  tend = 2*rpo(irpo).T;
  [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
  figure(1); set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend+10],'xtick',[0 .5 1]*tend); grid on;
  legend('b_1','c_1','b_2','c_2','b_3','c_3');
 axes(hax(2)); plot(tt,sum(aa.^2),'.-'); set(gca,'xlim',[0 tend+10],'xtick',[0 .5 1]*tend); grid on;
  
  %[tt, aa] = ksfmetd2(a0, L, h/10, tend, 1);
  %hold on; plot(tt,aa(1:6,:),'-');
%%
   ap = aa;  % Work in complex space
   for ia = 2:2, %size(aa,2), 
     va = ap(1:2:end,ia-1) + 1i*ap(2:2:end,ia-1);
     t = 1i*(1:N/2-1)'.*va;  P = (t*t')./(t'*t);
     da = ap(1:2:end,ia) + 1i*ap(2:2:end,ia) - va;
     dat = P*da;  dan = da - dat;  del = (t'*da)./(t'*t);
   end
%   

%%  
  ap = aa; % Work in real space
  G = zeros(N-2);  for ii = 1:2:N-2, G(ii:ii+1,ii:ii+1) = (ii+1)/2*[0 1;-1 0]; end
  %t = G*ap(:,1);  t = t./norm(t);  P = t*t';
  for ia = 2:2, %size(aa,2), 
    t = G*ap(:,ia-1);  P = (t*t')./(t'*t);
    da = ap(:,ia) - ap(:,ia-1);
    dat = P*da;  dan = da - dat;  del = (t'*da)./(t'*t);
    %phi = norm(dat)/2;
    %ap(:,ia-1) = ap(:,ia-1) + dat;
    %ap(:,ia:end) = ksfmshift2(ap(:,ia:end), -phi);
  end
  %hold on;  plot(tt, ap(1:6,:),'-');

%%
tic;
  G = zeros(N-2);  for ii = 1:2:N-2, G(ii:ii+1,ii:ii+1) = (ii+1)/2*[0 1;-1 0]; end
  for ii = 1:100000, t = G*ap(:,1); end
toc;

%%
tic;
  G = repmat(1:15,2,1); G(1,:) = -G(1,:);
  for ii = 1:100000, t1 = flipud(reshape(ap(:,1),2,[]).*G); t1 = t1(:); end
toc;

%% 
%   clear; load ks22f90h25t100; h = 0.25; ip = 7;
%   a0 = rpo(ip).a;  tend = rpo(ip).T;  na = length(a0);
%   [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%   figure(2); set(gcf,'pos',[350 650 700 300]); clf; 
%   plot(tt,aa(1:6,:),'o'); set(gca,'pos',[0.05 0.15 0.92 0.75],'xlim',[0 tend+5],'ygrid','on');
%   legend('b_1','c_1','b_2','c_2','b_3','c_3');
%   
%   G = zeros(na);  for ii = 1:2:na, G(ii:ii+1,ii:ii+1) = (ii+1)/2*[0 1;-1 0]; end
%   hh = 0.001;  nh = ceil(tend/hh);  
%   ah = zeros(length(a0),nh); ah(:,1) = a0; as = ah;
%   t = G*a0;  t = t./norm(t);  P = t*t';
%   for ih = 1:nh-1,
%     va = ksfm(0, ah(:,ih), L);  % KS flow 
%     ah(:,ih+1) = ah(:,ih) + hh.*va;
%    % t = G*ah(:,ih);  t = t./norm(t);  P = t*t';
%     vs = va - P*va;             % KS flow reduced by SO(2)
%     as(:,ih+1) = as(:,ih) + hh.*vs;
%   end
%   hold on;  plot((1:nh)'.*hh, ah(1:6,:), '-');  plot((1:nh)'.*hh, as(1:6,:), '--');
%   

%% Using Euler's method to integrate the flows
  clear; load ks22f90h25t100; h = 0.25;
  irpo = 5;  a0 = rpo(irpo).a; tend = 2*rpo(irpo).T; na = length(a0);
  h = 0.001;  nstp = ceil(tend/h);
  aa = zeros(na,nstp);  aa(:,1) = a0;  as = aa;  a4 = aa;
  G = zeros(na);  for ii = 1:2:na, G(ii:ii+1,ii:ii+1) = (ii+1)/2*[0 1;-1 0]; end  
  t = G*a0;  t = t./norm(t);  P = t*t';
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
  for istp = 1:nstp-1,
    va = ksfm(0, aa(:,istp), L);  % KS flow v(a)
    aa(:,istp+1) = aa(:,istp) + h*va;
   % a4(:,istp+1) = rk4(a4(:,istp), h, @(a)ksfm(0, a, L));
    t = G*as(:,istp);  t = t./norm(t);  P = t*t';
    vs = va - P*va;               % modified KS flow
    as(:,istp+1) = as(:,istp) + h*vs;
    if 0 && mod(istp-1,1000) == 0,
      [xx,uu] = ksfm2real(aa(:,istp), L, 64);  [xx,us] = ksfm2real(as(:,istp), L, 64);
      axes(hax(2));  plot(xx,[uu us],'.-'); grid on; pause;
    end
  end
 axes(hax(1));   plot((1:nstp)'*h,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(['RPO: $T = $~' num2str(tend/2) '~~$\ell = $~' num2str(rpo(irpo).s)],'interp','latex'); grid on;
  hold on; plot((1:nstp)'*h,a4(1:6,:),'o-'); 
%%  Using RK4 to integrate the flows
  clear; load ks22f90h25t100; h = 0.25;
  irpo = 5;  a0 = rpo(irpo).a; tend = 2*rpo(irpo).T; na = length(a0);
  h = 0.001;  nstp = ceil(tend/h);
  aa = zeros(na,nstp);  aa(:,1) = a0;  as = aa;  nta = zeros(nstp,1); nts = nta;
  G = zeros(na);  for ii = 1:2:na, G(ii:ii+1,ii:ii+1) = (ii+1)/2*[0 1;-1 0]; end  
  t = G*a0;  t = t./norm(t);  P = t*t';
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
  for istp = 1:nstp-1,
    va = ksfm(0, aa(:,istp), L);  % KS flow v(a)
    aa(:,istp+1) = rk4(aa(:,istp), h, @(a)ksfm(0, a, L));
    ta = G*aa(:,istp);  nta(istp+1) = ta'*(aa(:,istp+1)-aa(:,istp))./(ta'*ta);
    as(:,istp+1) = rk4(as(:,istp), h, @(a)ksfmperp(a, L, G));
    ts = G*as(:,istp);  nts(istp+1) = ts'*(as(:,istp+1)-as(:,istp))./(ts'*ts);
    if 0 && mod(istp-1,1000) == 0,
      [xx,uu] = ksfm2real(aa(:,istp), L, 64);  [xx,us] = ksfm2real(as(:,istp), L, 64);
      axes(hax(2));  plot(xx,[uu us],'.-'); grid on; pause;
    end
  end
 axes(hax(1));   plot((1:nstp)'*h,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(['RPO: $T = $~' num2str(tend/2) '~~$\ell = $~' num2str(rpo(irpo).s)],'interp','latex'); grid on;
  hold on; plot((1:nstp)'*h,as(1:6,:),'o-'); 
axes(hax(2));  plot((1:nstp)'*h, sum(aa.^2)-sum(as.^2), '.-'); grid on;
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
axes(hax(3));  plot((1:nstp)'*h, [nta nts],'.-');  grid on;
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
%%
% The RPO is shown over two periods.  In M/SO(2) it looks like a periodic
% orbit.  The shift variable correctly traces the shift of the RPO, 
% and so the difference between its values at the beginning and at the end of 
% the period is equal to the RPO shift (with a minus sign).
%%
  ippo = 6;  a0 = ppo(ippo).a1; tend = 2*ppo(ippo).T1;
  [tt,aa] = ksfmetd2(a0, L, h, tend, 1);
  set(gcf,'pos',[350 250 700 700]); clf; 
  hax = subplots(3,1,[0 0 0 0.05],[0.06 0.04 0.05 0.0]);
 axes(hax(1)); plot(tt,aa(1:6,:),'.-'); set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2);
  title(sprintf('PPO: $T = $%8.4f',tend/2),'interp','latex'); grid on;
  [af, ell] = ksfm_so2(aa, L);
 axes(hax(2)); plot(tt,ell*L,'.-'); 
  set(gca,'xlim',[0 tend],'ylim',[-L L]/2,'xtick',(0:2)*tend/2); grid on;
 axes(hax(3)); plot(tt,af(1:6,:),'.-'); 
  set(gca,'xlim',[0 tend],'xtick',(0:2)*tend/2); grid on;
%%
% In the case of a PPO, it is clear that it also becomes a periodic orbit in
% M/SO(2) and the shift variable shows the 'winding' and 'unwinding' of the
% shift over two periods.

%%
% Counter-example: why 'slicing' will not work for KS.
% Consider dynamical system defined on a space of periodic functions 
% in [-1, 1].  As with KS, we'll work with Fourier representation, with
% each Fourier mode given in polar coordinates: 
% (r_1, \theta_1, r_2, \theta_2, \ldots).

  a = zeros(30,1); a(1) = 1;  a(3) = 0.5;
  [x,u] = ksfm2real(a, 2); 
  figure(1); clf; plot(x,u,'.-'); hold on; grid on; p
  plot(x,ur,'r.-');
%  v = ifft(u([33:end 2:32]));
  
%% Counter-example: why 'slicing' will not work for KS.
  r1 = 1; r2 = 2;
  a0 = [r1; r2; 0]; t = 0;
  [t, a] = ode45(@(t,a)sliceflow(t,a,r1,r2), [0 1], a0);
  figure(1); clf; set(gcf,'paperpos',[1 10 19 6.5]); wysiwyg;
  hax = subplots(1,3,[0.02 0 0.16 0.12],[0.07 0.01 0 0]); 
  axes(hax(1)); plot(real(a(:,1)),imag(a(:,1)),'-',real(a(1,1)),imag(a(1,1)),'ro',real(a(end,1)),imag(a(end,1)),'k*'); grid on;
    xlabel('Re a_1'); ylabel('Im a_1');
  axes(hax(2)); plot(real(a(:,2)),imag(a(:,2)),'-',real(a(1,2)),imag(a(1,2)),'ro',real(a(end,2)),imag(a(end,2)),'k*'); grid on;
    title(['r_1 = ' num2str(r1) ',   r_2 = ' num2str(r2)], 'fontsize', 14);
    xlabel('Re a_2'); ylabel('Im a_2');
  axes(hax(3)); plot(t, a(:,3)./pi, '-'); grid on;
    xlabel('t'); ylabel('\theta/\pi');
  print -depsc2 sliceflow7.eps
  
%% Mapping unstable manifold of E2 to M/SO(2) using 1st Fourier mode
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  delta = 0.5*ere*period;
  a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
  [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
  [as, ss] = ksfm_so2(aa, L);
  figure(1); set(gcf,'pos',[782  89 620 855]); clf;
  for ii = 1:8, subplot(4,2,ii); 
    plot(tt,aa(2*ii,:),'.', tt,as(2*ii,:),'r.'); grid on; end;

%% Mapping unstable manifold of E2 to M/SO(2) using max_theta sum Im a_k
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 200;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  delta = 0.5*ere*period;
  a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
  [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
  kmax = 3; phase = (-0.98:0.001:1);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);  as = aa;  ss = zeros(1,size(aa,2));
  fac = exp(1i*pi*(1:size(v,1))'*phase);
  figure(1); clf;
  for ii = 1:1:size(aa,2),
    vs = fac.*repmat(v(:,ii),1,size(fac,2));
%    plot(phase,sum(imag(vs(1:3,:))),'.'); pause;
%    [vm, im] = max(imag(vs(1,:)));
    [vm, im] = max(sum(imag(vs(1:10,:))));
    as(1:2:end,ii) = real(vs(:,im)); as(2:2:end,ii) = imag(vs(:,im));
    ss(ii) = phase(im);
  end
  plot(tt,ss,'.');
  
%  [as, ss] = ksfm_so2(aa, L);
%  figure(1); set(gcf,'pos',[782  89 620 855]); clf;
%  for ii = 1:8, subplot(4,2,ii); 
%    plot(tt,aa(2*ii,:),'.', tt,as(2*ii,:),'r.'); grid on; end;

%% Tracking zeros and extremas of u(x,t)
  clear;  load kse22orbits;  h = 0.1;  tend = 150;
  switch 2,
    case 1,  a0 = rpo(37).a; 
    case 2,  k = 2; ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
      v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
      delta = 0.08*ere*period; a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);  
  end
  [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
  [x, uu] = ksfm2real(aa, L);
  figure(1); clf; set(gcf,'pos',[36   65   426   880]); 
  pcolor(x,tt,uu'); caxis([-3 3]); shading flat; hold on;
  ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
  xlabel('$x$','fontsize',14,'interp','latex');
  for ii = 1:size(aa,2),
    xz = ksfm2zeros(aa(:,ii), L);  plot(xz,tt(ii),'r.');
    xz1 = ksfm2zeros(aa(:,ii), L, 1);  plot(xz1,tt(ii),'w.');
    um = ksfm2anyx(aa(:,ii), L, xz1);
    [umi,imi] = min(um); [uma, ima] = max(um);
    plot(xz1(imi),tt(ii),'c.'); plot(xz1(ima),tt(ii),'m.'); 
    if uma > -umi, plot(xz1(ima),tt(ii),'k.');
    else plot(xz1(imi),tt(ii),'k.'); end
%    disp([xz1(imi) umi xz1(ima) uma]); pause;
  end

%% Mapping KSE solutions to M/O(2) using max/min or zeros of u(x)
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  delta = 0.5*ere*period;
  a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);  [x, ua] = ksfm2real(eq(k).a, L); 
  [tt, aa] = ksfmetd2(a0, L, h, tend, 1); 
  [x, uu, ux, uxx] = ksfm2real(aa, L);
  figure(1); clf; set(gcf,'pos',[36   265   426   680]); 
  pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
  ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
  xlabel('$x$','fontsize',14,'interp','latex');
  [ix,it] = find(uu(1:end-1,:)>0 & uu(2:end,:)<=0);  
  xz = zeros(size(ix));  uxz = xz;
  for ii = 1:length(ix),  
     ud = uu(ix(ii),it(ii))./(uu(ix(ii),it(ii))-uu(ix(ii)+1,it(ii)));
     xz(ii) = x(ix(ii)) + (x(ix(ii)+1)-x(ix(ii))).*ud;
     uxz(ii) = ux(ix(ii),it(ii)) + (ux(ix(ii)+1,it(ii))-ux(ix(ii),it(ii))).*ud;
  end
  
  
  hold on; plot(xz, tt(it), 'k.'); 
  
%% Plot coherent states and their eigenvectors in real space
 clear;  load kse22orbits;
 figure(1); clf; set(gcf,'pos',[489  33 695 916]); ne = 10;
 hax = subplots(ne/2+1,2,[0.01 0 0.02 0.03],[.04 .03 .03 .01]);
 switch 1, 
   case 1, % Equilibria
     k = 2;  a = eq(k).a;  eval = eq(k).eig(1:ne); evec = eq(k).evec(:,1:ne);
     axes(hax(1)); [xx,uu] = ksfm2real(a,L); plot(xx,uu,'.-'); grid on; axis tight;
     for ie = 1:ne,
       axes(hax(ie+2)); 
       if isreal(eval(ie)), re = evec(:,ie); 
         [xx,uu] = ksfm2real(re,L);  plot(xx,uu,'.-'); grid on; axis tight;
       else re = gsorth([real(evec(:,ie)) imag(evec(:,ie))]);
         [xx,uu] = ksfm2real(re(:,1),L);  plot(xx,uu,'.-'); grid on; hold on; 
         [xx,uu] = ksfm2real(re(:,2),L);  plot(xx,uu,'r.-'); grid on; axis tight;
       end, pause;
     end
   case 2, % RPOs
     k = 21; a = rpo(k).a; 
 end
 

