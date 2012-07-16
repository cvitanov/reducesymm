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
tp=104;
%Used for the projections
%Go to spectral space
%Layer1
ap=rv1(:,:,tp)';
ap=dst(ap);
ap=fft(ap')';
ap=ap(:,1:nx/2+1);
% Layer2
bp=rv2(:,:,tp)';
bp=dst(bp);
bp=fft(bp')';
bp=bp(:,1:nx/2+1);
%Magnitude of tangent vector
%Check Borders
%---------------------------%
aux_s=0;
for n=1:2
    if n==1
        up=ap;
    elseif n==2
        up=bp;
    end
    for k=1:nx/2
        k1=k+1;
        for l=0:ny-1
            l1=l+1;
            aux_s=aux_s+k^2*conj(up(l1,k1))*up(l1,k1);               
        end
    end
end
condition_s=real(-4*pi^2*aux_s/nx^2);

%Get a basis to project full state trajectories. In this case we project
%into a 2d space given by the template and the reflection about the x axis.
%Get e1
%---------------------------%
ea1=ap;
eb1=bp;
e01=[reshape(ea1,ny*(nx/2+1),1);reshape(eb1,ny*(nx/2+1),1)];
%Get e2
%---------------------------%
% Layer1
ap2=rv1(:,:,tp)';
ap2=-flipud(ap2);
ap2=dst(ap2);
ap2=fft(ap2')';
ap2=ap2(:,1:nx/2+1);
ea2=ap2;
% Layer2
bp2=rv2(:,:,tp)';
bp2=-flipud(bp2);
bp2=dst(bp2);
bp2=fft(bp2')';
bp2=bp2(:,1:nx/2+1);
eb2=bp2;
e02=[reshape(ea2,ny*(nx/2+1),1);reshape(eb2,ny*(nx/2+1),1)];
%Get the real projection vecotrs
e1=(e01+e02)/norm(e01+e02);
e2=(e01-e02)/norm(e01-e02);

%-------------------------------------------------------------------------%

%Interesting stuff wonder far away from the random inital conditions.
%-------------------------------------------------------------------------%
t1=0;
%Dimension of Arrays
a_rot=zeros(ny,nx/2+1);
b_rot=zeros(ny,nx/2+1);
rv1_rot=zeros(nx,ny,nt-tp);
rv2_rot=zeros(nx,ny,nt-tp);
rotf=zeros(nt-tp,1);
nextguess=zeros(nt-tp,1);
nextguess(t1+1)=0;
test=zeros(nt-tp,1);
condition_c=zeros(nt-tp,1);
condition_p=zeros(nt-tp,1);
condition=zeros(nt-tp,1);
ac=zeros(size(e1));
p1=zeros(nt-tp,1);
p2=zeros(nt-tp,1);
for t =tp:nt
    t1=t1+1;
    %Go to spectral space
    % Layer1
    a=rv1(:,:,t)';
    a=dst(a);
    a=fft(a')';
    a=a(:,1:nx/2+1);
    % Layer2
    b=rv2(:,:,t)';
    b=dst(b);
    b=fft(b')'; 
    b=b(:,1:nx/2+1);
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
    while stop==0 && e>0.5
        %Define the function value.
        [F,Fp]=sliceFunction2(a,ap,b,bp,rf,nx,ny);
        rf1=rf-F/Fp;
        e=abs(rf1-rf);
        %display(rf1)
        %If the rotaion frame change signs do not allow it.
        %         if rf1<0
        %             rf1=-rf1;
        %         end
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
%              display('A kick was given!')
%              display(rf)
%              display(rf1)
            rf=rf+15+5*rand();
            count=0;
            
        end
    end
    %A guess for the next rtf
    nextguess(t1+1)=rf+10*rand();
    rotf(t1)=rf;
    test(t1)=e;
    %Rotate to Slice
    %---------------------------%
    for k=0:nx/2
            k1=k+1;
        for l=0:ny-1
        l1=l+1;
            a_rot(l1,k1)=a(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
            b_rot(l1,k1)=b(l1,k1)*exp(1)^(2*pi*1i*k*rotf(t1)/nx);
        end
    end
    
    %Project to e1 and e2
    %---------------------------%
    ac=[reshape(a_rot,ny*(nx/2+1),1);reshape(b_rot,ny*(nx/2+1),1)];
    p1(t1)=real(dot(e1,ac));
    p2(t1)=real(dot(e2,ac));
    %Save Slices
    %---------------------------%
     a_rot=ifft(a_rot',nx,'symmetric');
     a_rot=idst(a_rot');
     b_rot=ifft(b_rot',nx,'symmetric');
     b_rot=idst(b_rot');
     rv1_rot(:,:,t1)=a_rot';
     rv2_rot(:,:,t1)=b_rot';
     
    %Check Borders
    %---------------------------%
    aux_c=0;
    aux_p=0;
    aux=0;
    for n=1:2
        if n==1
            u=a;
            up=ap;
        elseif n==2
            u=b;
            up=bp;
        end
        for k=1:nx/2
            k1=k+1;
            for l=0:ny-1
                l1=l+1;
                aux_c=aux_c+k^2*conj(u(l1,k1))*up(l1,k1)*exp(-2*pi*1i*rf*k/nx);
                aux_p=aux_p+k^2*conj(u(l1,k1))*u(l1,k1);

            end
        end
    end
    aux_c=aux_c+(nx/2)^2*conj(u(l1,k1))*up(l1,k1)*exp(-2*pi*1i*rf*k/nx);
    aux_p=aux_p+(nx/2)^2*k^2*conj(u(l1,k1))*u(l1,k1);
    
    condition_c(t1)=real(-4*pi^2*aux_c/nx^2);
    condition_p(t1)=real(-4*pi^2*aux_p/nx^2);
    condition(t1)=condition_c(t1)/(sqrt(condition_s)*sqrt(condition_p(t1)));
    %Run Info
    %---------------------------%   
    aux_txt=[num2str((t-tp)/(nt-tp)*100) '% Completed.' 'RotFrame:' num2str(rf1)];
    display(aux_txt);  
    clear a_rot b_rot ac
end
