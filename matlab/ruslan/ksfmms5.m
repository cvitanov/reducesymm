function [gm, dgm] = ksfmms5(am, L, h, mstp, np)
% function [gm, dgm] = ksfmms5(am, L, h, mstp, np);
%   5 - additional constraints: 
%       PFLG =  1: a(1) = 0, d a(2)/dt = 0 
%       PFLG = -1: PHASE = 0, a(1)*da(1)/dt + a(2)*da(2)/dt = 0

  global TT A F FT NFEVAL PHASE DD HH MSTP NP PFLG
 
  no = 2; scale = 1e-1;
  if nargin < 5, np = 1; end
  DD = L;  HH = h;  MSTP = mstp;  NP = np;
  na = round((length(am)-2)/np);  gm = zeros(size(am));
  tend = am(end-1);  PHASE = am(end);
  TT = [];  A = [];
  iaa = (1:na)+(np-1)*na; iar = (1:2:na)+(np-1)*na; iai = (2:2:na)+(np-1)*na;
  if nargout > 1, dtg = zeros(na,1); dgm = zeros(length(am)); end
  for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
    if nargout == 1,  [tt, aa] = ksfmstp2(am(ia), L, h, mstp, no);
    else [tt, aa, dgm(ia,ia)] = ksfmstp2(am(ia), L, h, mstp, no);
      dgm(ia,ia+na) = -eye(na); end
    gm(ia) = aa(:,end) - am((1:na)+ip*na);
    TT = [TT tt+h*(ip-1)*mstp]; A = [A aa]; end
  trem = tend-h*mstp*(np-1);
  if trem <= 0, disp('***Negative remaining integration time***'); pause; end
  if nargout == 1, [tt, aa] = ksfmetd2(am(iaa), L, h, trem, no);
  else [tt, aa, dgm(iaa,iaa), dtg] = ksfmetd2(am(iaa), L, h, trem, no); end
  TT = [TT tt+h*mstp*(np-1)]; A = [A aa];

  k = (-2i*pi/L).*(1:na/2)'; ek = exp(k.*PHASE);
  cg = ek.*(PFLG.*aa(1:2:end,end)+1i*aa(2:2:end,end));
  gm(iar) = real(cg) - am(1:2:na);
  gm(iai) = imag(cg) - am(2:2:na);
  
  cg = k.*cg; ddg = zeros(na,1);
  ddg(1:2:end-1) = real(cg);  ddg(2:2:end) = imag(cg);      % ddg
  
  gm(end-1:end) = 0.0;  
  
  [fa, dfa] = ksfm(0, am(1:na), L);
  if PFLG == 1,
    gm(end-1) = scale*am(1);  gm(end) = scale*fa(2);
  elseif PFLG == -1,
    gm(end-1) = scale*PHASE;  gm(end) = scale*(am(1).*fa(1) + am(2).*fa(2));
  end
    
  if nargout > 1,    
    cdg = repmat(ek,1,na).*(PFLG.*dgm(iar,iaa)+1i*dgm(iai,iaa));
    dgm(iar,iaa) = real(cdg); dgm(iai,iaa) = imag(cdg);
    dgm(iaa,1:na) = dgm(iaa,1:na) - eye(na);                % dag
    
    cf = ek.*(PFLG.*dtg(1:2:na)+1i*dtg(2:2:na));
    dtg(1:2:na) = real(cf); dtg(2:2:na) = imag(cf);
    dgm(iaa,end-1) = dtg;                                   % dtg
    dgm(iaa,end) = ddg;                                     % ddg
    dgm(end-1:end,:) = 0.0;
    
    if PFLG == 1,
      dgm(end-1,1) = scale*1.0;
      dgm(end,1:na) = scale*dfa(2,1:na);
    elseif PFLG == -1,
      dgm(end-1,end) = scale*1.0;
      dgm(end,1:na) = scale*(am(1).*dfa(1,1:na) + am(2).*dfa(2,1:na));
      dgm(end,1:2) = dgm(end,1:2) + scale*fa(1:2)';
    end
  end
  