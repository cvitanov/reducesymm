%Simulations Constants
% clc; clear;
cd('./run4')
nt=200;
nx=512*2;
ny=nx/4;
% load vort1.dat
% load vort2.dat
% rv1=reshape(vort1,nx,ny,nt);
% rv2=reshape(vort2,nx,ny,nt);
% save data_run9
% load data_run9

%% The real Thing. Reducing this symmetry.
%NOTE: NOT FUNCTIONAL AT THE MOMENT!!! ISSUES WITH THE ROTATING FRAME!!
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
nextguess(t1+1)=-100;
test=zeros(nt-tp,1);
condition=zeros(nt-tp,1);
for t =(tp+1):nt
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
    while e>0.000001
        aux1=0;
        aux2=0;
        aux=0;
        %Construct ap+Tap
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
                    aux=k*conj(u(l1,k1))*up(l1,k1)*exp(-2*pi*1i*rf*k/nx);
                    aux1=aux1+aux;
                    aux2=aux2+k*aux;
                end
            end
        end
        
        F1=-(2*pi*1i/nx)*aux1;
        F=real(F1);
        Fp=real(-(4*pi^2/nx^2)*aux2);
        rf1=rf-F/Fp;
        e=abs(rf1-rf);
        rf=rf1;
    end
    rotf(t1)=rf1;
    nextguess(t1+1)=rotf(t1);
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
                aux=aux+u(l1,k1)*up(l1,k1);
            end
        end
    end
    condition(t1)=-4*pi*aux/nx^2;
    %Run Info
    %---------------------------%   
    aux_txt=[num2str((t-tp)/(nt-tp)*100) '% Completed.'];
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