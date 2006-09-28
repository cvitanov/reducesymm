%% Exploring the property of full KSE with L = 22

%% Kuramoto-Sivashinsky Equation (KSE)
% $$ u_t = -uu_x - u_{xx} - u_{xxxx},  \quad x \in [-L/2, L/2]$$
%
% with periodic boundary condition: 
%
% $$ u(x+L,t) = u(x,t). $$
% 
% The equation has a reflection symmetry:
%
% $$ u(-x,t) = -u(x,t). $$

%% Properties of KSE equilibria with L = 22
% The KSE has a trivial equilibrium u = 0.  It is stable for L < 2\pi.  For
% L = 22 the eigenvalues are as follows:
  [h, dh] = ksfm(0, zeros(30,1), 22.0); [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e0eig = edh(ie);  e0vec = vdh(:,ie);
  disp(sprintf('%10.5f   (%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d, ... )\n',...
       [e0eig(1:10)'; e0vec(1:10,1:10)]));
%%
% Next to each eigenvalue above, I also list the corresponding eigenvector
% in the Fourier space:
% 
% $$ (\mathrm{Re}\tilde{u}_1, \mathrm{Im}\tilde{u}_1, \mathrm{Re}\tilde{u}_2, 
%     \mathrm{Im}\tilde{u}_2, \mathrm{Re}\tilde{u}_3, \mathrm{Im}\tilde{u}_3, \ldots) $$
%
% There are three pairs of unstable eigenvalues corresponding to the first
% three unstable modes: k = 2, 3, and 1, respectively.  For each mode, the 
% corresponding eigenvectors in the Fourier space lie in the plane 
%
% $$(\mathrm{Re} \tilde{u}_k, \mathrm{Im} \tilde{u}_k)$$
%
% Below are the KSE orbits starting close to u = 0 within each plane:
% 
% $$\tilde{u}_k(0) = (10^{-3}, 0)$$
% 
% Any other initial condition with the same magnitude within the plane 
% will yield the same orbit translated in x by d, where 
%
% $$ d = \frac{L}{2\pi k}\arg(\tilde{u}_k(0)) $$
  N = 32;  L = 22;  h = 0.1;  tend = 100;  np = 5;
  aa = cell(3,1);  uu = cell(3,1);  ux = cell(3,1);
  for k = [2 3 1],  a0 = zeros(2*N-2,1);  a0(2*k-1) = 1e-3;
    [tt, aa{k}] = ksfmedt(L, tend, a0, h, np);
    [x, uu{k}, ux{k}] = ksfm2real(aa{k}, L); end
  figure(1); clf; set(gcf,'pos',[5 450 700 400]);
  ax1 = axes('pos',[0.07 0.10 0.27 0.80]); pcolor(x,tt,uu{2}'); 
    shading flat; caxis([-3 3]); title('k = 2'); xlabel('x'); ylabel('t','rotat',0);
  ax2 = axes('pos',[0.39 0.10 0.27 0.80]); pcolor(x,tt,uu{3}'); 
    shading flat; caxis([-3 3]); title('k = 3'); xlabel('x');
  ax3 = axes('pos',[0.71 0.10 0.27 0.80]); pcolor(x,tt,uu{1}'); 
    shading flat; caxis([-3 3]); title('k = 1'); xlabel('x');
%%
% The k=2 and 3 perturbations converge to 2- and 3-wave equilibria (see
% below), respectively, while the k=1 perturbation also converges to the
% 2-wave, but shifted by L/8.
  load ks22uss1a;  [x, u1] = ksfm2real(ksfmshift(a0), d, 64);
  load ks22uss2a;  [x, u2] = ksfm2real(ksfmshift(a0), d, 64);
  load ks22uss3a;  [x, u3] = ksfm2real(ksfmshift(a0), d, 64);
  figure(1); clf; set(1,'position',[100 600 800 300]);
  ax1 = axes('pos',[0.06 0.12 0.28 0.74]);  plot(x,u1);  set(gca,'xlim',x([1 end]));
  title('1-wave equilibrium');  xlabel('x'); ylabel('u','rotat',0);
  ax2 = axes('pos',[0.38 0.12 0.28 0.74]);  plot(x,u2);  set(gca,'xlim',x([1 end]));
  title('2-wave equilibrium');  xlabel('x');
  ax3 = axes('pos',[0.70 0.12 0.28 0.74]);  plot(x,u3);  set(gca,'xlim',x([1 end]));
  title('3-wave equilibrium');  xlabel('x');
%%
% In addition to the trivial equilibrium, there are (at least ??) three
% equilibria corresponding to 1, 2, and 3-wave solutions (i.e. k-th Fourier 
% mode with the largest magnitude, k = 1, 2, 3).   The three equilibria are
% all symmetric with respect to the reflection symmetry.
  
%% Stability of equilibria
  load ks22uss1a;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e1eig = edh(ie);  e1vec = vdh(:,ie);
  load ks22uss2a;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e2eig = edh(ie);  e2vec = vdh(:,ie);
  load ks22uss3b;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e3eig = edh(ie);  e3vec = vdh(:,ie);
  disp(' 1-wave equilibrium  2-wave equilibrium  3-wave equilibrium');
  disp([e1eig(1:10) e2eig(1:10) e3eig(1:10)]);

%%
% The 1-wave equilibrium has two unstable planes.  The second plane is 
% in the odd subspace - the one found by Lan.  The 2- and 3-wave equilibria
% have only one unstable plane.  Therefore, it is expected that the 1-wave 
% equilibrium will have much less influence on the system evolution
% compared to the 2- and 3-wave equilibria. The 3-wave equilibrium has a 
% pair of real unstable eigenvalues which are equal to each other.  
% In this respect, it is similar to the eigenvalues of u = 0 equilibrium.

%% Evolution of orbits starting within the unstable planes of the equilibria
% Start with initial condition displaced by 0.001 from the equilibrium in
% the direction within the unstable plane(s). In the figures below the
% direction (along eigenvectors v_i) is indicated in the title.
  load ks22uss1a;  tend = 150;  np = 3;
  a1 = a0 + 1e-3.*real(e1vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  figure(1); set(gcf,'pos',[100 500 800 400]); clf;
  ax1 = axes('pos',[0.06 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 1: Re v_1');
  a1 = a0 + 1e-3.*real(e1vec(:,3));
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax2 = axes('pos',[0.30 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 2: Re v_3');
  load ks22uss2a;  tend = 150;  np = 2;
  a1 = a0 + 1e-3.*imag(e2vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax3 = axes('pos',[0.54 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('2-wave equilibrium: Im v_1');
  load ks22uss3b;  tend = 150;  np = 2;
  a1 = a0 + 1e-3.*e3vec(:,1);
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax4 = axes('pos',[0.78 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('3-wave equilibrium: v_1');
%%
  load ks22uss1a;  tend = 150;  np = 3;
  a1 = a0 - 1e-3.*real(e1vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  figure(1); set(gcf,'pos',[100 500 800 400]); clf;
  ax1 = axes('pos',[0.06 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 1: -Re v_1');
  a1 = a0 - 1e-3.*real(e1vec(:,3));
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax2 = axes('pos',[0.30 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 2: -Re v_3');
  load ks22uss2a;  tend = 150;  np = 2;
  a1 = a0 - 1.08e-3.*real(e2vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax3 = axes('pos',[0.54 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('2-wave equilibrium: -Re v_1');
  load ks22uss3b;  tend = 150;  np = 2;
  a1 = a0 + 1e-3.*e3vec(:,2);
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax4 = axes('pos',[0.78 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('3-wave equilibrium: v_2');
%%
% The 1-wave equilibrium in the 1st unstable plane falls straight into the
% chaotic attractor without passing close to any other equilibrium.  In the
% 2nd unstable plane, it evolves into a 2-wave equilibrium while staying in
% the odd subspace.  For a different initial displacement within the 2nd 
% plane, the orbit passes through a 3-wave equilibrium on the way to the 
% 2-wave equilibrium.  Thus it appears that the 2nd unstable plane (manifold) 
% of the 1-wave equilibrium intersects stable manifolds of both 2- and
% 3-wave equilibria.
%
% The 2-wave equilibrium falls into itself, but shifted by L/4.  
% It may also pass through the 3-wave equilibrium on the way to
% the L/4-shifted 2-wave equilibrium.
%
% The 3-wave equilibrium also evolves into the 2-wave equilibrium along 
% one of the unstable eigenvectors but falls into the
% chaotic attractor much faster (i.e. passes farther away from the 2-wave 
% equilibrium) along the other unstable eigenvector.
% 
% *Note:* The 2-wave and 3-wave plots above show *two* different mechanisms
% by which the 3-wave equilibrium evolves into the 2-wave one (see comments 
% below the 3-wave unstable manifold plot).

%% Investigating unstable manifolds of equilibria
% The 1- and 2-wave equilibria have complex pairs of unstable eigenvalues:
%
% $$ \lambda = \sigma \pm \mathrm{i}\omega, \quad \sigma > 0. $$
%
% An orbit starting close to the equilibrium will spiral out with a period
%
% $$ T = 2\pi/\omega. $$
%
% In order to trace out the unstable manifold, we start with a set of
% initial conditions
%
% $$\tilde{u} = \tilde{u}_{2w} + 10^{-4}\mathrm{e}^\delta v_1, \quad
%   \delta \in [0, \sigma T],$$
%
% where v_1 is a unit vector parallel to the real part of the unstable
% eigenvector.
%
  clear;  load kse22orbits;  k = 1;  p = 1;  h = 0.1;  tend = 75;
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,6))]);
  for delta = [0:0.05:ere*period 0.3 1.1],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend-delta/ere, a0, h, 2); 
    if delta == 0, av = v'*aa;
    else aa = [aa repmat(NaN,size(aa,1),size(av,2)-size(aa,2))]; av = [av; v'*aa]; end
    if delta == 0.3, tt1 = tt; aa1 = aa(:,1:length(tt1)); end
  end,
  figure(1); set(gcf,'pos',[100 400 800 500]); clf;
  ax1 = axes('pos',[0.35 0.06 0.62 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','k-');
    hold on; grid on; axis equal; view(-30,20);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color',[0 .8 0]);
    xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0);
  ax2 = axes('pos',[0.08 0.55 0.19 0.38]);
    [x, uu] = ksfm2real(aa1, L); 
    pcolor(x,tt1,uu'); caxis([-3 3]);
    shading flat;  ylabel('t','rotat',0);
    ht = title('Red orbit');  set(ht,'color','r','fontsize',15);
  ax3 = axes('pos',[0.08 0.08 0.19 0.38]);
    [x, uu] = ksfm2real(aa(:,1:length(tt)), L);  pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  ylabel('t','rotat',0); xlabel('x');
    ht = title('Green orbit');  set(ht,'color',[0 .8 0],'fontsize',15);
%%
% This is the manifold starting in the 1st unstable plane of the 1-wave
% equilibrium.  It falls directly into the chaotic attractor and does not 
% appear to contain any direct heteroclinic connections to other equilibria.
%
% Coordinate axes (v_1,v_2,v_3) are constructed by Gram-Schmidt
% orhonormalization of (Re e_1, Im e_1, Re e_6), where e_j are the
% eigenvectors of the 1-wave equilibrium.

%%
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 220;  av = [];
  ere = real(eq(k).eig(3));  period = 2*pi/imag(eq(k).eig(3));
  v = gsorth([real(eq(k).evec(:,3)) imag(eq(k).evec(:,3)) real(eq(k).evec(:,6))]);
  for delta = [0:0.05:ere*period 0.8 0.098 0.1028],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 2); av = [av; v'*aa];
    if delta == 0.8, aa1 = aa; end
    if delta == 0.098, aa2 = aa; end
  end,
  figure(1); set(gcf,'pos',[100 200 800 700]); clf;
  ax1 = axes('pos',[0.2 0.46 0.63 0.52]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','k-');
    hold on; grid on; axis equal; view(-115,15);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','b.-');
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color',[0 .8 0]);
    xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0);
  ax2 = axes('pos',[0.10 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa1, L); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  xlabel('x');  ylabel('t','rotat',0);
    ht = title('Blue orbit');  set(ht,'color','b','fontsize',15);
  ax3 = axes('pos',[0.40 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa2, L); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  xlabel('x');  ylabel('t','rotat',0);
    ht = title('Red orbit');  set(ht,'color','r','fontsize',15);
  ax4 = axes('pos',[0.70 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa, L);  pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  xlabel('x');  ylabel('t','rotat',0);
    ht = title('Green orbit');  set(ht,'color',[0 .8 0],'fontsize',15);
%%
% The manifold starting in the 2nd unstable plane of the 1-wave equilibrium
% converges to the 2-wave equilibrium.  It also contains a heteroclinic
% orbit to the 3-wave equilibrium.  Red and green orbits start close to
% this heteroclinic and then follow two different heteroclinics from 3- to
% 2-wave equilibrium (see 3-wave unstable manifold below).
%
% Coordinate axes (v_1,v_2,v_3) are constructed by Gram-Schmidt
% orhonormalization of (Re e_3, Im e_3, Re e_6), where e_j are the
% eigenvectors of the 1-wave equilibrium.

  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  for delta = [0:0.1:ere*period 2.7 0.24 0.254],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 2); av = [av; v'*aa];
    if delta == 2.7, aa1 = aa; end
    if delta == 0.24, aa2 = aa; end
  end,
  figure(1); set(gcf,'pos',[100 200 800 700]); clf;
  ax1 = axes('pos',[0.2 0.46 0.63 0.52]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','k-');
    hold on; grid on; axis equal; view(10,20);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','b.-');
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color',[0 .8 0]);
    xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0);
  ax2 = axes('pos',[0.10 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa1, L); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; xlabel('x'); ylabel('t','rotat',0);
    ht = title('Blue orbit');  set(ht,'color','b','fontsize',15);
  ax3 = axes('pos',[0.40 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa2, L); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; xlabel('x'); ylabel('t','rotat',0);
    ht = title('Red orbit');  set(ht,'color','r','fontsize',15);
  ax4 = axes('pos',[0.70 0.06 0.20 0.30]);
    [x, uu] = ksfm2real(aa, L);  pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; xlabel('x'); ylabel('t','rotat',0);
    ht = title('Green orbit');  set(ht,'color',[0 .8 0],'fontsize',15);
%%
% The unstable manifold of the 2-wave equilibrium converges to the 2-wave
% equilibrium shifted by L/4.  It also contains heteroclinic orbits to the
% 3-wave equilibrium.
%
% Coordinate axes (v_1,v_2,v_3) are constructed by Gram-Schmidt
% orhonormalization of (Re e_1, Im e_1, Re e_7), where e_j are the
% eigenvectors of the 2-wave equilibrium.

%%
% For the 3-wave equilibrium we have a pair of real unstable eigenvalues 
% equal to each other.  Therefore within the plane spanned by the 
% corresponding eigenvectors the orbits move radially away from the 
% equilibrium.  In order to trace out the unstable manifold, we start with 
% a set of initial conditions within the unstable plane
% 
% $$\tilde{u} = \tilde{u}_{3w} + 10^{-4}(v_1\cos\phi + v_2\sin\phi), \quad 
%   \phi \in [0, 2\pi],$$
%
% where v_1 and v_2 are orthonormal vectors within the plane spanned by the
% two unstable eigenvectors.  In the plot below, vector v_3 is parallel to 
% the eigenvector of the stable eigenvalue -0.4128.
  clear; load kse22orbits; k = 3; h = 0.1; delta = 1e-4; tend = 110; av = [];
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) eq(k).evec(:,4)]);
  for phi = [(0:0.1:6)*pi./3 1.3 3.395],
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmedt(L, tend, a0, h, 2);  av = [av; v'*aa];
    if phi == 1.3, aa1 = aa; end
  end,
  figure(1); set(gcf,'pos',[100 100 800 800]); clf;
  ax1 = axes('pos',[0.35 0.06 0.61 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','k-'); hold on; grid on;
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color',[0 .8 0]);
    xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0);
  ax2 = axes('pos',[0.08 0.69 0.19 0.25]);
    [x, uu] = ksfm2real(aa1, L); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  ylabel('t','rotat',0);
    ht = title('Red orbit');  set(ht,'color','r','fontsize',15);
  ax3 = axes('pos',[0.08 0.38 0.19 0.25]);
    [x, uu] = ksfm2real(aa, L);  pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  ylabel('t','rotat',0);
    ht = title('Green orbit');  set(ht,'color',[0 .8 0],'fontsize',15);
  ax4 = axes('pos',[0.08 0.07 0.19 0.25]);
    tend = 170;  phi = 0.25320172;
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmedt(L, tend, a0, h, 2);  av = v'*aa;
    [x, uu] = ksfm2real(aa, L);  pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat;  ylabel('t','rotat',0); xlabel('x');
    ht = title('Blue orbit');  set(ht,'color','b','fontsize',15);
  axes(ax1); plot3(av(1,:),av(2,:),av(3,:),'b.-'); axis equal; view(-150,30);
%%
% There are two types of heteroclinic connections from the 3- to the 
% 2-wave equilibrium.  They come out in the opposite directions from the
% 3-wave equilibrium and reach the 2-wave equlibrium with the same shift 
% (green and blue orbits).
% One of the heteroclinics (red and green orbits) is much more stable than
% the other one (blue).  But the blue one does occur in the heteroclinic
% networks 1,2-wave -> 3-wave -> L/4-shifted 2-wave shown above in the 1- and
% 2-wave unstable manifold plot.

%% Other representations for equilibria and orbits
% *8 Sep 2006, Predrag Cvitanovic wrote:*
% 
% |------------------------------------------------------------------|
% 
% A request - can you guys also plot all 3-dimensions for equilibria - 
% current eq (10) in rpo.tex, look at them by twirling them 
% 3-dimensionally? The $c$ integration constant also means something, 
% I am not quite sure what.
%
% These plots have advantage that they look the same for equilibria and
% relative equilibria (for relative equilibria point move in time, but the
% loop remains invariant). So we do not have the problem that from frame 
% of one of the equilibria rest look like circles - continuous symmetry is 
% built in.
%
% Perhaps there is some way of plotting the unstable manifolds too - I 
% would start with the hyperbolic ones rather than the complex planes - 
% that would enable us to look at heteroclinic connections. Fourier reps 
% are a bit wierd.
%
% |------------------------------------------------------------------|
  load kse22orbits; colr = ['b';'r';'k'];
  figure(1); set(gcf,'pos',[420 500 560 450]); clf;
  ax1 = axes('Position',[0.05 0.0828 0.788 0.907]);
  axis(ax1,[-2.597 2.597 -3.002 1.459 -2.804 2.804]);
  ylabel(ax1,'u_x'); zlabel(ax1,'u_{xx}','rotat',0); xlabel(ax1,'u');
  view(ax1,[25 25]); grid(ax1,'on'); hold(ax1,'all');
  for k = 1:3, [x,u,ux,uxx] = ksfm2real(eq(k).a,L,128);
    plot3(u,ux,uxx,['.-' colr(k)]); hold on; end
  grid on; axis equal; set(gca,'pos',[0.03 0.05 0.9 0.94]); view(25,25);
  set(ax1,'pos',[0.02 0.08 0.79 0.91]);  set(get(gca,'xlabel'),'pos',[16.6 -39 15]);  
  set(get(gca,'ylabel'),'pos',[20 -37 15.5]);  set(get(gca,'zlabel'),'pos',[-1.9 -5.8 1.6]);
  lg1 = legend(ax1,{'1-wave','2-wave','3-wave'});  set(lg1,'pos',[0.75 0.3 0.18 0.14]);  
%%
% Another way to view the equilibria is in the space
%
% $$ (u, u_x, u_{xx}) $$
%
% where they look like non-intersecting closed loops parameterized by x.
% 
%
%% 
% *Ruslan's reply:*
% I've plotted the equilibria in coordinates of Eq. (10).  Note that, as 
% x changes from 0 to L, the loops representing 2- and 3-wave 
% equilibria are travelled 2 and 3 times, respectively.  I
% doubt we will be able to get much mileage out of this representation.  
% Even though it looks nice, the problem is that these coordinates are not 
% giving a proper phase-space representation - phase space points look like 
% closed loops here.
%
% I'm getting more convinced now that we cannot get away from Fourier
% representation.  We just need to better understand how to exploit the
% fact that our phase space is actually composed of complex planes.
% One idea, which is similar to what you proposed, is to look at the
% Fourier modes in polar coordinates:
%
% $$ \tilde{u}_k = r_k\mathrm{e}^{-\mathrm{i}k\theta_k}$$
%
% It is convenient to use this form because translation
%
% $$ u(x) \rightarrow u(x+d) $$
%
% is represented by a phase shift
%
% $$ \theta_k \rightarrow \theta_k + d/\tilde{L} $$
%
% which is the same for all modes.  If we don't care about the phase
% information (i.e. if we don't want to see equilibria as circles) then we
% simply plot things in radial coordinates (e.g. r_1, r_2, r_3).  Below are
% some examples.
  load kse22orbits; colr = ['b';'r';'k']; delta = 1e-4; tend = 170;
  figure(1); set(gcf,'pos',[420 500 560 450]); clf;
  ax1 = axes('Position',[0.03 0.11 0.9 0.85]);
  view(ax1,[100 10]); grid(ax1,'on'); hold(ax1,'all');
  for k = 1:3, 
    r1 = abs(eq(k).a(1)+1i*eq(k).a(2));  r2 = abs(eq(k).a(3)+1i*eq(k).a(4));
    r3 = abs(eq(k).a(5)+1i*eq(k).a(6));  plot3(r1,r2,r3,['*' colr(k)]); hold on; end
  grid on; axis equal; axis([0 1 0 1 0 1.5]);
  xlabel(ax1,'r_1'); ylabel(ax1,'r_2'); zlabel(ax1,'r_3','rotat',0);
  set(get(gca,'xlabel'),'pos',[0.5 1.2 0]);  set(get(gca,'ylabel'),'pos',[1.6 0.5 0]);
  set(get(gca,'zlabel'),'pos',[1.2 -0.1 0.75]);
  v = gsorth([real(eq(3).evec(:,1)) imag(eq(3).evec(:,1))]);
  for phi = 0.25320172 + [pi 0],
    a0 = eq(3).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmedt(L, tend, a0, h, 1);
    r1 = abs(aa(1,:)+1i*aa(2,:));  r2 = abs(aa(3,:)+1i*aa(4,:));
    r3 = abs(aa(5,:)+1i*aa(6,:));
    plot3(r1,r2,r3,'.-','color',[0 .8 0]);
  end,
%%
% These are plots of the two heteroclinic orbits from E3 to E2.