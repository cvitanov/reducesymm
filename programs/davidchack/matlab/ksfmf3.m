function g = ksfmf3(a, d, h);
%  function g = ksfmf3(a, d, h)  -  Associated flow for the full KSFM
%  3 - for use with fsolve

global TT A F FT NFEVAL hp1 hp2

  NFEVAL = NFEVAL + 1;  if NFEVAL > 32, np = 2; else np = 0; end  

  [TT, A, F, FT] = ksfmedt(d, a(end), a(1:end-1), h, np);
  if np == 0,  g = A - a(1:end-1);  else  g = A(:,end) - a(1:end-1); end
  g = [g; F'*g];
  
  if NFEVAL > 32,
    figure(3);  set(hp1,'xdata',A(1,:),'ydata',A(2,:)); 
    plot(A(1,1),A(2,1),'ro',A(1,end),A(2,end),'ko');
    figure(5);  v = A(1:2:end,:) + 1i*A(2:2:end,:);
    v = [zeros(1,size(A,2)); v; zeros(1,size(A,2)); flipud(conj(v))];  
    u = real(fft(v));  set(hp2,'ydata',TT,'zdata',u','cdata',u'); axis tight;
    NFEVAL = 0; end
return;
