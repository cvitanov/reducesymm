function as = ksfmshift(a)
%  Shift the phase of FMs to get real(u_k) = 0 for smallest k 
%  such that |u_k| > 0.

  v = a(1:2:end) + 1i*a(2:2:end);  mv = abs(v);
  iv = min(find(mv > 1e-6));
  avk = (angle(v(iv))-pi/2)./iv;
  vs = exp(-1i*avk.*(1:length(v))').*v;
  as = a;  as(1:2:end) = real(vs);  as(2:2:end) = imag(vs);
  