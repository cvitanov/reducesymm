%Simulations Constants
% clc; clear;
% cd('./run4')
% path(path,'../')
% nt=200;
% nx=512*2;
% ny=nx/4;
% load vort1.dat
% load vort2.dat
% rv1=reshape(vort1,nx,ny,nt);
% rv2=reshape(vort2,nx,ny,nt);
% save data_run4
% load data_run4
%% The real Thing. Reducing this symmetry.
%Select the template (ap)
%-------------------------------------------------------------------------%
tp=100;
%Go to spectral space
% Layer1
a1p=rv1(:,:,tp)';
a1p=dst(a1p);
a1p=fft(a1p')';
% Layer2
b1p=rv2(:,:,tp)';
b1p=dst(b1p);
b1p=fft(b1p')';
%Magnitude of tangent vector
%Check Borders
%---------------------------%
aux_s=0;
for n=1:2
    if n==1
        up=a1p;
    elseif n==2
        up=b1p;
    end
    for k=1:nx-1
        k1=k+1;
        for l=0:ny-1
            l1=l+1;
            aux_s=aux_s+k^2*conj(up(l1,k1))*up(l1,k1);
        end
    end
end
condition_s=real(-4*pi^2*aux_s/nx^2);

%-------------------------------------------------------------------------%

%Interesting stuff wonder far away from the random inital conditions.
%-------------------------------------------------------------------------%
t1=0;
%Dimension of Arrays
a1_rot=zeros(ny,nx);
b1_rot=zeros(ny,nx);
rv1_rot=zeros(nx,ny,nt-tp);
rv2_rot=zeros(nx,ny,nt-tp);
rotf=zeros(nt-tp,1);
nextguess=zeros(nt-tp,1);
nextguess(t1+1)=0;
test=zeros(nt-tp,1);
condition_c=zeros(nt-tp,1);
condition_p=zeros(nt-tp,1);
condition=zeros(nt-tp,1);
for t =tp:nt
    t1=t1+1;
    %Go to spectral space
    % Layer1
    a1=rv1(:,:,t)';
    a1=dst(a1);
    a1=fft(a1')';
    % Layer2
    b1=rv2(:,:,t)';
    b1=dst(b1);
    b1=fft(b1')';  
    %Start our Newton Raphson
    %---------------------------%
    %Chose a initial guess for the rotating frame
    rf=nextguess(t1);
    e=1; 
    stop=0;
    count=0;
    while stop==0
        %Define the funvtion value.
        F=sliceFunction(a1,a1p,b1,b1p,rf,nx,ny);
        Fp1=sliceFunction(a1,a1p,b1,b1p,rf+1,nx,ny);
        Fm1=sliceFunction(a1,a1p,b1,b1p,rf-1,nx,ny);
        Fp=(Fp1-Fm1)/2;
        %Next rotating frame
        rf1=rf-round(F/Fp);
        e=abs(rf1-rf);
        %Looping count
        if rf+1==rf1 || rf-1==rf1
            count=count+1;
        end
        %Stop condition
        if rf==rf1 || count>3;
            stop=1;
        end
        rf=rf1;        
    end
    nextguess(t1+1)=rf1;
    rotf(t1)=rf1;
    test(t1)=e;
    %Rotate to Slice
    %---------------------------%
    for k=0:nx-1
            k1=k+1;
        for l=0:ny-1
        l1=l+1;
            a1_rot(l1,k1)=a1(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
            b1_rot(l1,k1)=b1(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
        end
    end
    %Save Slices
     a1_rot=ifft(a1_rot');
     a1_rot=idst(a1_rot');
     b1_rot=ifft(b1_rot');
     b1_rot=idst(b1_rot');
     rv1_rot(:,:,t1)=a1_rot';
     rv2_rot(:,:,t1)=b1_rot';

    
    %Check Borders
    %---------------------------%
    aux_c=0;
    aux_p=0;
    aux=0;
    for n=1:2
        if n==1
            u=a1;
            up=a1p;
        elseif n==2
            u=b1;
            up=b1p;
        end
        for k=1:nx-1
            k1=k+1;
            for l=0:ny-1
                l1=l+1;
                aux_c=aux_c+k^2*conj(u(l1,k1))*up(l1,k1);
                aux_p=aux_p+k^2*conj(u(l1,k1))*u(l1,k1);
            end
        end
    end
    condition_c(t1)=real(-4*pi^2*aux_c/nx^2);
    condition_p(t1)=real(-4*pi^2*aux_p/nx^2);
    condition(t1)=condition_c(t1)/(condition_s*condition_p(t1));
    %Run Info
    %---------------------------%   
    aux_txt=[num2str((t-tp)/(nt-tp)*100) '% Completed.' 'RotFrame:' num2str(rf1)];
    display(aux_txt);
    
end

%Deleted Code
%-------------------------------------------------------------------------%
%N1
%Construct ap+Tap. This term is only imaginary.
% aux=0;
% for n=1:2
%     if n==1
%         u=a1p;
%     elseif n==2
%         u=b1p;
%     end
%     for k=0:nx-1
%         k1=k+1;
%         for l=0:ny-1
%             l1=l+1;
%             aux=aux+k*abs(u(l1,k1))^2;
%         end
%     end
% end
%
% F2=-aux*2*pi*i/nx;

%N2
% %Test in 2D
% for t=115:1:nt
%
%  A=rv1(:,:,t)';
%  B=dst(A);
%  C=fft(B')';
%
%  b1=rv2(:,:,t)';
%  b1=dst(b1);
%  b1=fft(b1')';
%
%  %Shifting the image
% for n=0:ny-1
%     n1=n+1;
%     for m=0:nx-1
%     m1=m+1;
%     C2(n1,m1)=C(n1,m1)*exp(1)^(-2*pi*i*m*200/nx);
%     end
% end
%  subplot(2,1,1)
%  imagesc(A)
%
%  subplot(2,1,2)
%  imagesc(b1)
%
% %  subplot(2,2,3)
% %  D2=ifft(C2');
% %  E2=idst(D2');
% %  imagesc(real(E2))
% %
% %  subplot(2,2,4)
% %  D=ifft(C');
% %  E=idst(D');
% %  imagesc(real(E))
% end

%N3
% %% Test in 1D
% clc; clear;
% a=zeros(16,1);
% a(4:8)=1;
% b=fft(a);
% subplot(2,2,1)
% plot(a)
% axis([0 16 0 1.5])
% subplot(2,2,2)
% plot(abs(b))
% %Shifting the image
% for n=0:15
%     n1=n+1;
%     b2(n1,1)=b(n1)*exp(1)^(-2*pi*i*n*3/16);
% end
% b2(1)=real(b2(1));
% b2(9)=real(b2(9));
% a1=ifft(b2);
% subplot(2,2,3)
% plot(real(a1))
% axis([0 16 0 1.5])
%
% subplot(2,2,4)
% a=ifft(b)
% plot(a)
% axis([0 16 0 1.5])
