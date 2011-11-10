function gi = ksfmfns(ai, d, h, ni);
% function gi = ksfmfns(ai, d, h, ni);

global TT A F FT NFEVAL PHASE DD

  nn = round((length(ai)-2)/(ni(1)+1));  gi = ai;  
  TT = [];  A = [];  PHASE = ai(end);  DD = d;
  for ii = 1:ni(1),
    [tt, aa] = ksfmedt(d, h*(ni(2)+.5), ai((1:nn)+(ii-1)*nn), h, 1);
    TT = [TT tt(1:end-1)+h*(ii-1)*ni(2)]; A = [A aa(:,1:end-1)];
    gi((1:nn)+(ii-1)*nn) = aa(:,end-1) - ai((1:nn)+ii*nn);  
  end
  [tt, aa, F, FT] = ksfmedt(d, ai(end-1), ai((1:nn)+ni(1)*nn), h, 1);
  TT = [TT tt+h*ni(1)*ni(2)];  A = [A aa];
  
  ek = exp(-2i*pi/d.*PHASE.*(1:nn/2)');
  v = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
  
  gi((1:2:nn)+ni(1)*nn) = real(v) - ai(1:2:nn);
  gi((2:2:nn)+ni(1)*nn) = imag(v) - ai(2:2:nn);

  vf = ek.*(FT(1:2:end-1)+1i*FT(2:2:end));  ftp = zeros(size(FT));  
  ftp(1:2:end-1) = real(vf);  ftp(2:2:end) = imag(vf);
  
  v = 1i*(1:nn/2)'.*v;  app = zeros(nn,1);
  app(1:2:end-1) = real(v);  app(2:2:end) = imag(v);

  gi(end-1:end) = [ftp'; app']*gi((1:nn)+ni(1)*nn);
  
return;