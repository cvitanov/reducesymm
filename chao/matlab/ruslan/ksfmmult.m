function [gm, dgm] = ksfmmult(am, L, nstp, np);
% function [gm, dgm] = ksfmmult(am, L, nstp, np);
  
%  clear; load kse22orbits; ipo = 10; np = 1;
%  nstp = round(rpo(ipo).T./0.25);  h = rpo(ipo).T./nstp;
%  mstp = round(nstp/np); am = rpo(ipo).a(1:30);  aa = rpo(ipo).a(1:30);
%  for ip = 1:np-1, [tt, aa] = ksfmedt(L, mstp*h, aa, h); am = [am; aa]; end
%  am = [am; h; rpo(ipo).d];  nargin = 4;  nargout = 2;

  global TT A F FT NFEVAL PHASE DD
%  persistent noutput
 
  alpht = 1e-2;  alphd = 1e-3; no = 2;
  if nargin < 4, np = 1; end
  na = round((length(am)-2)/np);  gm = zeros(size(am));
  mstp = round(nstp/np); h = am(end-1)/nstp; PHASE = am(end); DD = L;
  TT = []; A = [];  dtg = zeros(na*np,1);
  if nargout > 1, ddtg = zeros(na*np, na); dgm = zeros(length(am)); end
  for ip = 1:np-1, ia = (1:na)+(ip-1)*na;  rstp = mstp./nstp;
    if nargout == 1,  [tt, aa] = ksfmstp(am(ia), L, h, mstp, no);
      dtg(ia) = ksfm(0, aa(:,end), L);  dtg(ia) = dtg(ia).*rstp;
    else, [tt, aa, dgm(ia,ia)] = ksfmstp(am(ia), L, h, mstp, no);
      [dtg(ia), ddtg(ia,:)] = ksfm(0, aa(:,end), L);  dtg(ia) = dtg(ia).*rstp;
      dgm(ia,ia+na) = -eye(na);  ddtg(ia,:) = ddtg(ia,:).*rstp; end
    gm(ia) = aa(:,end) - am((1:na)+ip*na);
    TT = [TT tt+h*(ip-1)*mstp]; A = [A aa]; end
  ia = (1:na)+(np-1)*na;  rstp = (nstp-mstp*(np-1))./nstp;
  if nargout == 1,  [tt, aa] = ksfmstp(am(ia), L, h, nstp-mstp*(np-1), no);
    dtg(ia) = ksfm(0, aa(:,end), L);  dtg(ia) = dtg(ia).*rstp;
  else, [tt, aa, dgm(ia,ia)] = ksfmstp(am(ia), L, h, nstp-mstp*(np-1), no); 
    [dtg(ia), ddtg(ia,:)] = ksfm(0, aa(:,end), L);  dtg(ia) = dtg(ia).*rstp;
    ddtg(ia,:) = ddtg(ia,:).*rstp; end
  TT = [TT tt+h*(np-1)*mstp];  A = [A aa];

  iar = (1:2:na)+(np-1)*na; iai = (2:2:na)+(np-1)*na;
  k = (-2i*pi/L).*(1:na/2)'; ek = exp(k.*PHASE);
  cg = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
  gm(iar) = real(cg) - am(1:2:na);  gm(iai) = imag(cg) - am(2:2:na);

  cf = ek.*(dtg(iar)+1i*dtg(iai));  fnp = dtg(ia);
  dtg(iar) = real(cf); dtg(iai) = imag(cf);                 % dtg
  
  cg = k.*cg; ddg = zeros(na,1);
  ddg(1:2:end-1) = real(cg);  ddg(2:2:end) = imag(cg);      % ddg
  
  gm(end-1) = -alpht.*(dtg'*gm(1:end-2));
  gm(end) = -alphd.*(ddg'*gm(ia));
 
  if nargout > 1,
    datg = zeros(np*na,na);
    for ip = 1:np, ia = (1:na)+(ip-1)*na;
      datg(ia,:) = ddtg(ia,:)*dgm(ia,ia); end
    ia = (1:na)+(np-1)*na; 
    cdg = repmat(ek,1,na).*(datg(iar,:)+1i*datg(iai,:));
    datg(iar,:) = real(cdg);  datg(iai,:) = imag(cdg);      % datg
    
    cdg = repmat(ek,1,na).*(dgm(iar,ia)+1i*dgm(iai,ia));
    dgm(iar,ia) = real(cdg); dgm(iai,ia) = imag(cdg);
    dgm(ia,1:na) = dgm(ia,1:na) - eye(na);                  % dag
        
    cdg = repmat(k,1,na).*cdg;  dadg = zeros(na);
    dadg(1:2:na,:) = real(cdg); dadg(2:2:na,:) = imag(cdg); % dadg

    d2tg = zeros(np*na,1);
    for ip = 1:np-1, ia = (1:na)+(ip-1)*na;
      d2tg(ia) = ddtg(ia,:)*dtg(ia); end
    ia = (1:na)+(np-1)*na;  d2tg(ia) = ddtg(ia,:)*fnp;
    cdtg = ek.*(d2tg(iar)+1i*d2tg(iai));
    d2tg(iar) = real(cdtg);  d2tg(iai) = imag(cdtg);        % d2tg
    
    cf = k.*cf;  dtdg = zeros(na,1);
    dtdg(1:2:end) = real(cf);  dtdg(2:2:end) = imag(cf);    % dtdg
    
    cg = k.*cg;  d2dg = zeros(na,1);
    d2dg(1:2:end) = real(cg);  d2dg(2:2:end) = imag(cg);    % d2dg
    
    dgm(1:na*np,end-1) = dtg;                               % dtg
    dgm(ia,end) = ddg;                                      % ddg
    
    dgm(end-1,1:na*np) = -alpht.*(dtg'*dgm(1:end-2,1:end-2));
    for ip = 1:np, ia = (1:na)+(ip-1)*na;
      dgm(end-1,ia) = dgm(end-1,ia) - ...
        alpht.*(gm(ia)'*datg(ia,:)); end                    % dagt
    
    ia = (1:na)+(np-1)*na;
    dgm(end,ia) = -alphd.*(ddg'*dgm(ia,ia) + gm(ia)'*dadg);
    if np > 1, dgm(end,1:na) = dgm(end,1:na) + alphd.*ddg'; end   % dagd
    
    dgm(end-1,end-1) = -alpht.*(dtg'*dtg + gm(1:na*np)'*d2tg);
    dgm(end,end-1) = -alphd.*(ddg'*dtg(ia) + gm(ia)'*dtdg);
    dgm(end-1,end) = -alpht.*(ddg'*dtg(ia) + gm(ia)'*dtdg);
    dgm(end,end) = -alphd.*(ddg'*ddg + gm(ia)'*d2dg);
  end
%  NFEVAL = NFEVAL + 1; if NFEVAL == 1; noutput = 0; end
%  if noutput < NFEVAL, disp(num2str([norm(dtg) norm(ddg) am(end-1:end)']));
%    noutput = noutput + 100; end