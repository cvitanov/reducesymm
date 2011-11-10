function adot = sliceflow(t, a, r1, r2)
%  Example of simple linear flow in on 2 Fourier modes reduced to the slice
%  Note that a(1:2) are complex, while a(3) = \theta

  om1 = pi/2;  om2 = 3*pi;

  adot = zeros(3,1); 
  adot(3) = (om1*r1*2*real(a(1)) + 2*om2*r2*2*real(a(2)))./(r1*2*real(a(1)) + 4*r2*2*real(a(2)));
  adot(1) = 1i*om1*a(1) - adot(3)*1i*a(1);
  adot(2) = 1i*om2*a(2) - adot(3)*2i*a(2);
  