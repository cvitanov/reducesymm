function [gm, dgm] = ksfmmult2(am, L, nstp, np);
% function [gm, dgm] = ksfmmult2(am, L, nstp, np);
  
%  clear; load kse22orbits; ipo = 10; np = 1;
%  nstp = round(rpo(ipo).T./0.25);  h = rpo(ipo).T./nstp;
%  mstp = round(nstp/np); am = rpo(ipo).a(1:30);  aa = rpo(ipo).a(1:30);
%  for ip = 1:np-1, [tt, aa] = ksfmedt(L, mstp*h, aa, h); am = [am; aa]; end
%  am = [am; h; rpo(ipo).d];  nargin = 4;  nargout = 2;

  global TT A F FT NFEVAL PHASE DD
 
  alpht = 1e-0;  alphd = 1e-0; no = 2;  dh = 1e-7;
  if nargin < 4, np = 1; end
  na = round((length(am)-2)/np);  gm = zeros(size(am));
  mstp = round(nstp/np); h = am(end-1)/nstp; PHASE = am(end); DD = L;
  dh = h*dh; TT = []; A = []; f = zeros(na*np,1);
  iaa = (1:na)+(np-1)*na; iar = (1:2:na)+(np-1)*na; iai = (2:2:na)+(np-1)*na;
  if nargout > 1, df = zeros(na*np, na); dtg = zeros(na*np,1);
    dgm = zeros(length(am)); end
  for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
    if nargout == 1,  [tt, aa] = ksfmstp(am(ia), L, h, mstp, no);
      f(ia) = ksfm(0, am(ia+na), L);
    else, [tt, aa, dgm(ia,ia)] = ksfmstp(am(ia), L, h, mstp, no);
      [f(ia), df(ia,:)] = ksfm(0, am(ia+na), L);
      [ttd, dtg(ia)] = ksfmstp(am(ia), L, h+dh, mstp, 0);
      dtg(ia) = (dtg(ia)-aa(:,end))./(nstp*dh); dgm(ia,ia+na) = -eye(na); end
    gm(ia) = aa(:,end) - am((1:na)+ip*na);
    TT = [TT tt+h*(ip-1)*mstp]; A = [A aa]; end
  if nargout == 1,  [tt, aa] = ksfmstp(am(iaa), L, h, nstp-mstp*(np-1), no);
    f(iaa) = ksfm(0, am(1:na), L);
  else, [tt, aa, dgm(iaa,iaa)] = ksfmstp(am(iaa), L, h, nstp-mstp*(np-1), no); 
    [f(iaa), df(iaa,:)] = ksfm(0, am(1:na), L); 
    [ttd, dtg(iaa)] = ksfmstp(am(iaa), L, h+dh, nstp-mstp*(np-1), 0);
    dtg(iaa) = (dtg(iaa)-aa(:,end))./(nstp*dh); end
  TT = [TT tt+h*(np-1)*mstp]; A = [A aa];

  k = (-2i*pi/L).*(1:na/2)'; ek = exp(k.*PHASE);
  cg = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
  gm(iar) = real(cg) - am(1:2:na);  gm(iai) = imag(cg) - am(2:2:na);
  
  cg = k.*cg; ddg = zeros(na,1);
  ddg(1:2:end-1) = real(cg);  ddg(2:2:end) = imag(cg);      % ddg
  
  gm(end-1) = -alpht.*(f'*gm(1:end-2));
  gm(end) = -alphd.*(ddg'*gm(iaa));
 
  if nargout > 1,
%    datg = zeros(np*na,na);
%    for ip = 1:np, ia = (1:na)+(ip-1)*na;
%      datg(ia,:) = ddtg(ia,:)*dgm(ia,ia); end
%    ia = (1:na)+(np-1)*na; 
%    cdg = repmat(ek,1,na).*(datg(iar,:)+1i*datg(iai,:));
%    datg(iar,:) = real(cdg);  datg(iai,:) = imag(cdg);      % datg
    
    cdg = repmat(ek,1,na).*(dgm(iar,iaa)+1i*dgm(iai,iaa));
    dgm(iar,iaa) = real(cdg); dgm(iai,iaa) = imag(cdg);
    dgm(iaa,1:na) = dgm(iaa,1:na) - eye(na);                 % dag
        
    cdg = repmat(k,1,na).*cdg;  dadg = zeros(na);
    dadg(1:2:na,:) = real(cdg); dadg(2:2:na,:) = imag(cdg);  % dadg

%    d2tg = zeros(np*na,1);
%    for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
%      d2tg(ia) = ddtg(ia,:)*dtg(ia); end
%    ia = (1:na)+(np-1)*na;  d2tg(ia) = ddtg(ia,:)*fnp;
%    cdtg = ek.*(d2tg(iar)+1i*d2tg(iai));
%    d2tg(iar) = real(cdtg);  d2tg(iai) = imag(cdtg);        % d2tg
    
%    cf = k.*cf;  dtdg = zeros(na,1);
%    dtdg(1:2:end) = real(cf);  dtdg(2:2:end) = imag(cf);    % dtdg
    
    cg = k.*cg;  d2dg = zeros(na,1);
    d2dg(1:2:end) = real(cg);  d2dg(2:2:end) = imag(cg);    % d2dg
    
    cf = ek.*(dtg(iar)+1i*dtg(iai));
    dtg(iar) = real(cf); dtg(iai) = imag(cf);
    dgm(1:na*np,end-1) = dtg;                               % dtg
    
    dgm(iaa,end) = ddg;                                     % ddg

    dgm(end-1,1:na*np) = -alpht.*(f'*dgm(1:end-2,1:end-2));
    for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
      dgm(end-1,ia+na) = dgm(end-1,ia+na) -...
        alpht.*(gm(ia)'*df(ia,:)); end                      % dagt
    dgm(end-1,1:na) = dgm(end-1,1:na) - alpht.*(gm(iaa)'*df(iaa,:));
    
    dgm(end,iaa) = -alphd.*(ddg'*dgm(iaa,iaa) + gm(iaa)'*dadg);
    if np > 1, dgm(end,1:na) = dgm(end,1:na) + alphd.*ddg'; end % dagd
    
    ddtg = dtg(iaa);
    cf = k.*(ddtg(1:2:end)+1i*ddtg(2:2:end));
    ddtg(1:2:end) = real(cf);  ddtg(2:2:end) = imag(cf);
    
    dgm(end-1,end-1) = -alpht.*(f'*dtg);
    dgm(end,end-1) = -alphd.*(ddg'*dtg(iaa) + gm(iaa)'*ddtg);
    dgm(end-1,end) = -alpht.*(ddg'*f(iaa));
    dgm(end,end) = -alphd.*(ddg'*ddg + gm(iaa)'*d2dg);
  end
%  NFEVAL = NFEVAL + 1; if NFEVAL == 1; noutput = 0; end
%  if noutput < NFEVAL, disp(num2str([norm(dtg) norm(ddg) am(end-1:end)']));
%    noutput = noutput + 100; end