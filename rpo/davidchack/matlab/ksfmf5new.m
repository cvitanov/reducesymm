function [g, jac] = ksfmf5new(t, a, d, h, C, alpha);
%  function g = ksfmf5(t, a, d, h, C, alpha)  -  Associated flow for the full KSFM
%  5 - for quasiperiodic orbits

global TT A F FT NFEVAL PHASE DD
  
  na = length(a)-2;
  if nargout == 1, NFEVAL = NFEVAL + 1;
    [TT, A] = ksfmedt(d, a(end-1), a(1:end-2), h, 2);
    FT = ksfm(t, A(:,end), d);
  else,  NFEVAL = NFEVAL + na + 3;
    [TT, A, DA] = ksfmjedt(d, a(end-1), a(1:end-2), h, 2);
    [FT, DF] = ksfm(t, A(:,end), d);
  end
  F = ksfm(t, a(1:end-2), d);  PHASE = a(end);  DD = d;
%  [F, DF] = ksfm(t, a(1:end-2), d);  PHASE = a(end);  DD = d;
  
  k = (-2i*pi/d).*(1:na/2)';  ek = exp(k.*PHASE);
  cg = ek.*(A(1:2:end,end)+1i*A(2:2:end,end));
  
  g = a(1:end-2);  g(1:2:end) = real(cg) - g(1:2:end);
  g(2:2:end) = imag(cg) - g(2:2:end);
  
  cf = ek.*(FT(1:2:end-1)+1i*FT(2:2:end));  dtg = zeros(na,1);
  dtg(1:2:end) = real(cf);  dtg(2:2:end) = imag(cf);
  
  cg = k.*cg;  ddg = zeros(na,1);
  ddg(1:2:end) = real(cg);  ddg(2:2:end) = imag(cg);
  
  if nargout > 1,
    dug = eye(na);
    cdg = repmat(ek,1,na).*(DA(1:2:end-1,:)+1i*DA(2:2:end,:));
    dug(1:2:end,:) = real(cdg)-dug(1:2:end,:);
    dug(2:2:end,:) = imag(cdg)-dug(2:2:end,:);
    
    ddug = eye(na);
    cdg = repmat(k,1,na).*cdg;
    ddug(1:2:end,:) = real(cdg);  ddug(2:2:end,:) = imag(cdg);
    
    dutg = DF*DA;
    cdg = repmat(ek,1,na).*(dutg(1:2:end-1,:)+1i*dutg(2:2:end,:));
    dutg(1:2:end,:) = real(cdg);  dutg(2:2:end,:) = imag(cdg);
    
    dt2g = DF*FT;
    cf = ek.*(dt2g(1:2:end-1,:)+1i*dt2g(2:2:end,:));
    dt2g(1:2:end) = real(cf);  dt2g(2:2:end) = imag(cf);
    
    cf = k.*cf;  ddtg = zeros(na,1);
    ddtg(1:2:end) = real(cf);  ddtg(2:2:end) = imag(cf);
    
    cg = k.*cg;  dd2g = zeros(na,1);
    dd2g(1:2:end) = real(cg);  dd2g(2:2:end) = imag(cg);

    jac = [C; -alpha(1).*dtg'; -alpha(2).*ddg']*[dug dtg ddg];
    jac(na+1,:) = jac(na+1,:) - alpha(1).*g'*[dutg dt2g ddtg];
    jac(na+2,:) = jac(na+2,:) - alpha(2).*g'*[ddug ddtg dd2g];
%    jac = C*[dug dtg ddg];
%    jac = [jac; [1 zeros(1,na+1)]; 2.*[F(1).*DF(1,:)+F(2).*DF(2,:) 0 0]];
  end
  
  g = [C; -alpha(1).*dtg'; -alpha(2).*ddg']*g; 
%  g = [C*g; a(1); F(1).^2+F(2).^2];
return;
