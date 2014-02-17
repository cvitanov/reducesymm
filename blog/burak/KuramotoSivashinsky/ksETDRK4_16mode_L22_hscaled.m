%Adaptation of kursiv.m from Kassam ^ Trefethen (2005)

%Spatial grid:
d = 22;
N = 32;
%Initial condition:
%u = cos(4*pi*x/d).*(1+sin(4*pi*x/d));

rpo = [1.081769255879405645e-01 0.000000000000000000e+00  -1.130645021532460798e-01 ...  
	2.735234400271993951e-02  -2.300369007695817619e-02  2.743365567776075153e-02 ...
	4.242109469805525057e-01  -3.221375201839944691e-02  3.517350195620121411e-01 ...
	4.196323928512094570e-01  7.405822221667555938e-02  -4.911698645198029345e-01 ...
	-2.619084037151255262e-01  8.869647954573157966e-03  2.667068425090810685e-02 ...
	-1.245857190912121742e-01  1.848625450932936676e-01  -1.886910780372257068e-01 ...
	-4.364329592632312099e-02  -8.603322827952401136e-03  -4.893648586116418342e-02 ...
	-4.227361593906937137e-02  -5.743046880645331920e-02  6.141174799787345318e-02 ...
	3.556320072237982056e-03  -2.647610106987533310e-02  -3.295731184608969265e-03 ... 
	-1.760410218329051119e-02  -1.449156681256755577e-02  1.551927466950007474e-02]';
    
rpo = [3.70955584e-01   1.59178104e-16  -5.35414753e-01   6.72233738e-01 ...
  -1.60102241e-01   2.77492997e-01   2.26119482e-01   9.20202541e-02 ...
   1.01234154e-01  -7.61348435e-02   2.86999279e-03  -5.71633872e-02 ...
  -2.19911067e-02  -8.45432215e-03  -9.70953085e-03   5.74786687e-03 ...
   4.19544211e-04   4.07227658e-03   1.61695688e-03   7.51275352e-04 ...
   5.51797781e-04  -3.99001588e-04  -2.26058837e-05  -2.53441206e-04 ...
  -9.51845395e-05  -3.49971067e-05  -2.92711613e-05   2.18049965e-05 ...
   5.41318900e-07   1.12892181e-05]';


%v = zeros(32,1);
%v(2:16) = rpo(1:2:29) + 1i*rpo(2:2:30);
%v(32:-1:18) = rpo(1:2:29) - 1i*rpo(2:2:30);
%v(17)=0;

x0 = rpo;
v = [0; x0(1:2:end-1)+1i*x0(2:2:end); 0; x0(end-1:-2:1)-1i*x0(end:-2:2)];
u=real(ifft(v));

%Template:
slicep = zeros(N,1);
slicep(2) = 1;
slicep(N) = 1;

%ETDRK4 scalars:
h0 = 0.1;                        % time step
k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers
T = 1i * [0:N/2-1 0 -N/2+1:-1]'; %U(1) generator
%Linear term:  
Lold = k.^2 - k.^4;
%L = Lold + T;           
%L = Lold;
L = Lold;

E0 = exp(h0*L); E20 = exp(h0*L/2);
M = 16;                         % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);    % roots of unity
LR0 = h0*L(:,ones(M,1)) + r(ones(N,1), :);
Q0  = h0*real(mean(           (exp(LR0/2) - 1)./LR0                 ,2));
f10 = h0*real(mean(   (-4-LR0+exp(LR0).*(4-3*LR0+LR0.^2))./LR0.^3   ,2));
f20 = h0*real(mean(       (2+LR0+exp(LR0).*(-2+LR0))./LR0.^3        ,2));
f30 = h0*real(mean(   (-4-3*LR0-LR0.^2+exp(LR0).*(4-LR0))./LR0.^3   ,2));

%Time-stepping loop:
uu = u; tt = 0; vv = v;
tmax = 2*32.80617425;

%Nonlinear term:
g = 0.5i*k*N;
Nold = @(x) g.*fft(real(ifft(x)).^2); fft(real(ifft(v)).^2);

vold = @(x) Lold.*x + Nold(x);

%Reduced velocity:
tp = T.*slicep;
vhat = @(x) vold(x) - ((vold(x)'*tp)/((T.*x)'*tp)) * (T.*x);

%N = @(x) Nold(x) + (real(x(2))-1)*vold(x) - (imag(vold(x)(2)) + 1)*(T.*x);
%N = @(x) Nold(x) + (real(x(2))-1)*vold(x) - imag(vold(x)(2))*(T.*x);
N = @(x) Nold(x) - (imag(vold(x)(2))/real(x(2)))*(T.*x);

vnew = @(x) L.*x + N(x);

t = 0;

while t < tmax
    
    if abs(real(v(2))) < 1

        h = h0*abs(real(v(2)));
        
        E = exp(h*L); E2 = exp(h*L/2);
        
        LR = (h/h0)*LR0;
        Q  = h*real(mean(           (exp(LR/2) - 1)./LR              ,2));
        f1 = h*real(mean(   (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3   ,2));
        f2 = h*real(mean(       (2+LR+exp(LR).*(-2+LR))./LR.^3       ,2));
        f3 = h*real(mean(   (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3   ,2));
        
    else
        
        h = h0;
        E = E0; E2 = E20;

        LR = LR0;
        Q  = Q0;
        f1 = f10;
        f2 = f20;
        f3 = f30;
    
    end

    t = t+h;

    Nv = N(v);
    a = E2.*v + Q.*Nv;
    Na = N(a);
    b = E2.*v + Q.*Na;
    Nb = N(b);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = N(c);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    %if mod(n,nplt)==0
        u=real(ifft(v));
        uu = [uu, u]; tt = [tt, t]; vv = [vv, v];
    %end
end

figure(1)
plot3(real(vv(2,:)), real(vv(3,:)), imag(vv(3,:)))
xlabel('x_1')
ylabel('x_2')
zlabel('y_2')

#save('data/sspsolution.dat', 'ta')
box off

print -dpng ksslice_hscaled.png
