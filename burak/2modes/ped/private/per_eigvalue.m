
function pereig=per_eigvalue(M)

    n=size(M,1); m=size(M,2)/n;  pereig=[]; 
    j=1;   
    while j<n
        %fprintf(1,'j=%d\n',j);       
        if M(j+1,(m-1)*n+j)==0, mhess=0; sig=1; % deal with real eigenvalue.
            for k=1:m,
                nor=abs(M(j,(k-1)*n+j)); sig=sig*sign(M(j,(k-1)*n+j)); mhess=log(nor)+mhess;  
            end 
            pereig=[pereig;[mhess,sig]];
            j=j+1;
        else  % complex conjugate pair
            mhess=eye(2); mag=0;
            for k=1:m
                mhess=M(j:j+1,(k-1)*n+j:(k-1)*n+j+1)*mhess;
                nor=max(max(abs(mhess))); mag=mag+log(nor); mhess=mhess/nor;
            end
            %normalize the eigenvalue. move the norm to the exponent part.
            cplx=eig(mhess); mag2=abs(cplx); cplx=cplx./mag2; mag=mag+log(mag2);            
            pereig=[pereig; [mag,cplx]]; 
            j=j+2;
        end
    end
    % deal with the last diagonal element of Mm.
    if M(n,(m-1)*n+n-1)==0,   mhess=0; sig=1;
        for k=1:m
            nor=abs(M(n,k*n)); sig=sig*sign(M(n,(k*n))); mhess=log(nor)+mhess;  
        end 
        pereig=[pereig;[mhess,sig]];
    end

end


 
