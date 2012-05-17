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
a1p=a1p(:,1:nx/2+1);
% Layer2
b1p=rv2(:,:,tp)';
b1p=dst(b1p);
b1p=fft(b1p')';
b1p=b1p(:,1:nx/2+1);
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
    for k=1:nx/2
        k1=k+1;
        for l=0:ny-1
            l1=l+1;
            aux_s=aux_s+k^2*conj(up(l1,k1))*up(l1,k1);
        end
    end
end
condition_s=real(4*pi^2*aux_s/nx^2);

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
    a=rv1(:,:,t)';
    a=dst(a);
    a=fft(a')';
    a1=a(:,1:nx/2+1);
    % Layer2
    b=rv2(:,:,t)';
    b=dst(b);
    b=fft(b')'; 
    b1=b(:,1:nx/2+1);
    %Start our Newton Raphson
    %---------------------------%
    %Chose a initial guess for the rotating frame
    rf=nextguess(t1);
    e=2; 
    stop=0;
    count=0;
    count2=0;
    ea=0;
    sum=0;
    while stop==0 && e>0.1
        %Define the function value.
        [F,Fp]=sliceFunction2(a1,a1p,b1,b1p,rf,nx,ny);
        rf1=rf-F/Fp;
        e=abs(rf1-rf);
        %display(rf1)
        %If the rotaion frame change signs do not allow it. This is what
        %might be reffered to something hardwired to the code. As I know
        %that the first values of rf1 will be possitive.
        if rf1<0
            rf1=-rf1;
        end
        %Check if the solution is the minimum and trow it to a expected
        %interval if it is not there.
        if t1>1           
            if rf1<nextguess(t1-1)-100 || rf1>nextguess(t1-1)+100
                rf1=nextguess(t1-1)+15+5*rand();
            end
        end
        %If more than a rotation is made return to an equivalent value of
        %rf1
        if rf1>nx
        rf=rf1-floor(rf1/(nx))*nx;
        else
        rf=rf1;
        end
        %Give it a kick if the search is endless
        count=count+1;
        if count>20
            rf=rf+15+5*rand();
            count=0;
        end
    end
    %Backward finite difference scheme to find the next rtf
    nextguess(t1+1)=round(rf)+10;
    rotf(t1)=rf;
    test(t1)=e;
    %Rotate to Slice
    %---------------------------%
    for k=0:nx/2
            k1=k+1;
        for l=0:ny-1
        l1=l+1;
            a1_rot(l1,k1)=a1(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
            b1_rot(l1,k1)=b1(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
        end
    end
     
    %Save Slices
     a1_rot=ifft(a1_rot',nx,'symmetric');
     a1_rot=idst(a1_rot');
     b1_rot=ifft(b1_rot',nx,'symmetric');
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
        for k=1:nx/2
            k1=k+1;
            for l=0:ny/2
                l1=l+1;
                aux_c=aux_c+k^2*conj(u(l1,k1))*up(l1,k1);
                aux_p=aux_p+k^2*conj(u(l1,k1))*u(l1,k1);
            end
        end
    end
    condition_c(t1)=real(4*pi^2*aux_c/nx^2);
    condition_p(t1)=real(4*pi^2*aux_p/nx^2);
    condition(t1)=condition_c(t1)/(sqrt(condition_s)*sqrt(condition_p(t1)));
    %Run Info
    %---------------------------%   
    aux_txt=[num2str((t-tp)/(nt-tp)*100) '% Completed.' 'RotFrame:' num2str(rf1)];
    display(aux_txt);
    
end
