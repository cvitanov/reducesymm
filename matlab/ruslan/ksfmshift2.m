function as = ksfmshift2(a, phi)
%  Shift the phase of KS FMs by angle phi

  v = a(1:2:end,:) + 1i*a(2:2:end,:);
  vs = repmat(exp(-1i*phi.*(1:size(v,1))'),1,size(v,2)).*v;
  as = a;  as(1:2:end,:) = real(vs);  as(2:2:end,:) = imag(vs);
  