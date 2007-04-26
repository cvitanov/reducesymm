function [gm, dgm] = ksfmms(am, L, h, mstp, np);
% function [gm, dgm] = ksfmms(am, L, h, mstp, np);

  global TT A F FT NFEVAL PHASE DD HH MSTP NP PFLG
 
  alpht = 1e-4;  alphd = 1e-2;  no = 2;
  if nargin < 5, np = 1; end
  DD = L;  HH = h;  MSTP = mstp;  NP = np;
  na = round((length(am)-2)/np);  gm = zeros(size(am));
  tend = am(end-1);  PHASE = am(end);
  TT = [];  A = [];  f = zeros(na*np,1);
  iaa = (1:na)+(np-1)*na; iar = (1:2:na)+(np-1)*na; iai = (2:2:na)+(np-1)*na;
  if nargout > 1, df = zeros(na*np, na); dtg = zeros(na,1);
    dgm = zeros(length(am)); end
  for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
    if nargout == 1,  [tt, aa] = ksfmstp(am(ia), L, h, mstp, no);
      f(ia) = ksfm(0, am(ia+na), L);
    else, [tt, aa, dgm(ia,ia)] = ksfmstp(am(ia), L, h, mstp, no);
      [f(ia), df(ia,:)] = ksfm(0, am(ia+na), L);
      dgm(ia,ia+na) = -eye(na); end
    gm(ia) = aa(:,end) - am((1:na)+ip*na);
    TT = [TT tt+h*(ip-1)*mstp]; A = [A aa]; end
  trem = tend-h*mstp*(np-1);
  if trem <= 0, disp('***Negative remaining integration time***'); pause; end
  if nargout == 1,  [tt, aa] = ksfmetd(am(iaa), L, h, trem, no);
    f(iaa) = ksfm(0, am(1:na), L);
  else, [tt, aa, dgm(iaa,iaa), dtg] = ksfmetd(am(iaa), L, h, trem, no); 
    [f(iaa), df(iaa,:)] = ksfm(0, am(1:na), L); end
  TT = [TT tt+h*mstp*(np-1)]; A = [A aa];

  k = (-2i*pi/L).*(1:na/2)'; ek = exp(k.*PHASE);
  cg = ek.*(PFLG.*aa(1:2:end,end)+1i*aa(2:2:end,end));
  gm(iar) = real(cg) - am(1:2:na);  gm(iai) = imag(cg) - am(2:2:na);
  
  cg = k.*cg; ddg = zeros(na,1);
  ddg(1:2:end-1) = real(cg);  ddg(2:2:end) = imag(cg);    % ddg
  
  gm(end-1) = -alpht.*(f'*gm(1:end-2));
  gm(end) = -alphd.*(ddg'*gm(iaa));
 
  if nargout > 1,    
    cdg = repmat(ek,1,na).*(PFLG.*dgm(iar,iaa)+1i*dgm(iai,iaa));
    dgm(iar,iaa) = real(cdg); dgm(iai,iaa) = imag(cdg);
    dgm(iaa,1:na) = dgm(iaa,1:na) - eye(na);                 % dag
        
    cdg = repmat(k,1,na).*cdg;  dadg = zeros(na);
    dadg(1:2:na,:) = real(cdg); dadg(2:2:na,:) = imag(cdg);  % dadg
    
    cg = k.*cg;  d2dg = zeros(na,1);
    d2dg(1:2:end) = real(cg);  d2dg(2:2:end) = imag(cg);    % d2dg
    
    cf = ek.*(PFLG.*dtg(1:2:na)+1i*dtg(2:2:na));
    dtg(1:2:na) = real(cf); dtg(2:2:na) = imag(cf);
    dgm(iaa,end-1) = dtg;                                   % dtg
    
    dgm(iaa,end) = ddg;                                     % ddg
    
    dgm(end-1,1:na*np) = -alpht.*(f'*dgm(1:end-2,1:end-2));
    for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
      dgm(end-1,ia+na) = dgm(end-1,ia+na) -...
        alpht.*(gm(ia)'*df(ia,:)); end                      % dagt
    dgm(end-1,1:na) = dgm(end-1,1:na) - alpht.*(gm(iaa)'*df(iaa,:));
    
    dgm(end,iaa) = -alphd.*(ddg'*dgm(iaa,iaa) + gm(iaa)'*dadg);
    if np > 1, dgm(end,1:na) = dgm(end,1:na) + alphd.*ddg'; end % dagd
    
    ddtg = dtg;  cf = k.*(ddtg(1:2:end)+1i*ddtg(2:2:end));
    ddtg(1:2:end) = real(cf);  ddtg(2:2:end) = imag(cf);
    
    dgm(end-1,end-1) = -alpht.*(f(iaa)'*dtg);
    dgm(end,end-1) = -alphd.*(ddg'*dtg + gm(iaa)'*ddtg); %%% ???
    dgm(end-1,end) = -alpht.*(ddg'*f(iaa));
    dgm(end,end) = -alphd.*(ddg'*ddg + gm(iaa)'*d2dg);
  end
  