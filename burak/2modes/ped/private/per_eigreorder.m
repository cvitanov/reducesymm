
function [pereig,Ev]=per_eigreorder(M,Q)
%function [pereig,Ev]=xdper_eig(M,Q) calculates the eigenvalues and 
%eigenvectors of sequence of matrices M. 
%input:
%       M: periodic schur form of matrices. M1, M2,...Mm. the first m-1
%           matrices are upper-triangular. Mm is quasi-upper triangular with
%           a 1x1 block for each real eigenvalue and a 2x2 block for each pair
%           of conjugate eigenvalues.
%       Q: orthogonal matrices.
%output:
%       pereig: eigenvalues.
%               left column is the logrithm of magnitude.
%               right column is the sign of eigenvalue. For complex number,
%               it is normalized.
%       Ev: eigenvectors along the oribt.
    n=size(M,1); m=size(M,2)/n;  pereig=[]; Ev=[];
    j=1;   
    while j<n
        %fprintf(1,'j=%d\n',j);       
        if M(j+1,(m-1)*n+j)==0, mhess=0; sig=1; % deal with real eigenvalue.
            for k=1:m,
                nor=abs(M(j,(k-1)*n+j)); sig=sig*sign(M(j,(k-1)*n+j)); mhess=log(nor)+mhess;  
            end 
            pereig=[pereig;[mhess,sig]];
            %calculate real eigenvectors
            ve=zeros(n,m);
            v=eigv(M,j-1,1);
            for i=1:m, k=mod(i-2,m)+1; ve(:,i)=Q(:,(k-1)*n+1:k*n)*v(1:n,i); end 
            % index sequence: m,1,2,...,m-1
            Ev=[Ev;ve]; j=j+1;            
        else  % complex conjugate pair
            mhess=eye(2); mag=0;
            for k=1:m
                mhess=M(j:j+1,(k-1)*n+j:(k-1)*n+j+1)*mhess;
                nor=max(max(abs(mhess))); mag=mag+log(nor); mhess=mhess/nor;
            end
            %normalize the eigenvalue. move the norm to the exponent part.
            cplx=eig(mhess); mag2=abs(cplx); cplx=cplx./mag2; mag=mag+log(mag2);            
            pereig=[pereig; [mag,cplx]]; 
            %calculate complex eigenvectors
            v=eigv(M,j-1,2);
            ve=zeros(n,m);
            for i=1:m, k=mod(i-2,m)+1; ve(:,i)=Q(:,(k-1)*n+1:k*n)*v(1:n,i); end
            Ev=[Ev;ve]; 
            for i=1:m, k=mod(i-2,m)+1; ve(:,i)=Q(:,(k-1)*n+1:k*n)*v(n+1:2*n,i); end
            Ev=[Ev;ve];
            j=j+2;
        end
    end
    % deal with the last diagonal element of Mm.
    if M(n,(m-1)*n+n-1)==0,   mhess=0; sig=1;
        for k=1:m
            nor=abs(M(n,k*n)); sig=sig*sign(M(n,(k*n))); mhess=log(nor)+mhess;  
        end 
        pereig=[pereig;[mhess,sig]];
        %calculate real eigenvectors
        ve=zeros(n,m);v=eigv(M,n-1,1);
        for i=1:m, k=mod(i-2,m)+1; ve(:,i)=Q(:,(k-1)*n+1:k*n)*v(1:n,i); end
        Ev=[Ev;ve];
    end


    end


function v=eigv(M,p1,flag)
%function v=eigv(M,p1,flag)
%   gets the (p1+1)th eigenvector by reorder (1:P1,1:P1) block with (p1+1, p1+1)
%   for real case or (p1+1:p1+2, p1+1,p1+2) for complex case.
%input: 
%   M: seqence of (quasi-)uppper triangular matrices
%   p1: position
%   flag: flag=1 for real case. otherwise for complex case.
%output:
%       v:(p1+1)th eigenvector along the orbit.     
n=size(M,1); m=size(M,2)/n;
if flag==1, % real eigenvector.
    if p1==0, v=zeros(n,m); v(1,:)=1;
    else
        [pse,t12]=sp_pse(M,p1,flag);       
        spparms('autommd',0);
        x=pse\t12;
        spparms('autommd',1);
        v=zeros(n,m); v(1:p1,:)=reshape(x,p1,m); v(p1+1,:)=1; ...
          v=normc(v);    
    end
else % complex eigenvector pair
    mhess=eye(2); j=p1+1;
    for k=1:m,
        mhess=M(j:j+1,(k-1)*n+j:(k-1)*n+j+1)*mhess;
        nor=max(max(abs(mhess))); mhess=mhess/nor;
    end
    [ve,~]=eig(mhess); v1=zeros(n,m); v2=zeros(n,m);
    
    if p1~=0
        [pse,t12]=sp_pse(M,p1,flag);       
        spparms('autommd',0);
        x=pse\t12;
        spparms('autommd',1);
        x=reshape(x,p1,2*m);
    end
    
    for i=1:m
        vt=zeros(n,2); vt(p1+1:p1+2,:)=eye(2); if p1~=0, vt(1:p1,:)=x(:,(i-1)*2+1:i*2); end
        vt=vt*ve; v1(:,i)=vt(:,1)/norm(vt(:,1)); v2(:,i)=vt(:,2)/norm(vt(:,2));
        ve=M(j:j+1,(i-1)*n+j:(i-1)*n+j+1)*ve; % use rotated 2x2 eigenvector.
        ve(:,1)=ve(:,1)/norm(ve(:,1)); ve(:,2)=ve(:,2)/norm(ve(:,2));
    end
    v=[v1;v2];
    
end

end

function [pse,t12]=sp_pse(M,p1,flag)
% build the sparse matrix 'pse' and dense column vector 't12';

n=size(M,1); m=size(M,2)/n;
if flag==1,
    t12=zeros(p1*m,1); bs=p1^2+p1;
    ii=zeros(bs*m,1);jj=zeros(bs*m,1);ss=zeros(bs*m,1);
    for i=1:m
        t12((i-1)*p1+1:i*p1)=-M(1:p1,(i-1)*n+p1+1);
            
        xx=((i-1)*p1+1:i*p1)'; if i~=m, yy=(i*p1+1:(i+1)*p1)'; else yy=(1:p1)'; end,   
        it=kron(ones(p1,1),xx); jt=kron(xx,ones(p1,1));
        is1=reshape(M(1:p1,(i-1)*n+1:(i-1)*n+p1),p1^2,1); is2=-M(p1+1,(i-1)*n+p1+1)*ones(p1,1);
        ii((i-1)*bs+1:i*bs)=[it;xx];
        jj((i-1)*bs+1:i*bs)=[jt;yy]; 
        ss((i-1)*bs+1:i*bs)=[is1;is2]; 
    end
    pse=sparse(ii,jj,ss,p1*m,p1*m);
else
    t12=zeros(2*p1*m,1); bs=2*p1^2+4*p1;
    ii=zeros(bs*m,1);jj=zeros(bs*m,1);ss=zeros(bs*m,1);
    for i=1:m
         t12((i-1)*2*p1+1:i*2*p1)=-M(1:p1,(i-1)*n+p1+1:(i-1)*n+p1+2);
         
         xx=((i-1)*2*p1+1:i*2*p1)'; if i~=m, yy=(i*2*p1+1:(i+1)*2*p1)'; else yy=(1:2*p1)'; end,
         it1=kron(ones(p1,1),xx(1:p1)); it2=kron(ones(p1,1),xx(p1+1:2*p1));
         jt1=kron(xx(1:p1),ones(p1,1)); jt2=kron(xx(p1+1:2*p1),ones(p1,1));
         is1=reshape(M(1:p1,(i-1)*n+1:(i-1)*n+p1),p1^2,1);
         R22=-M(p1+1:p1+2,(i-1)*n+p1+1:(i-1)*n+p1+2)';
         is2=kron(R22(:),ones(p1,1)); 
         
         ii((i-1)*bs+1:i*bs)=[it1;it2;xx(1:p1);xx(p1+1:2*p1);xx(1:p1);xx(p1+1:2*p1)];
         jj((i-1)*bs+1:i*bs)=[jt1;jt2;yy(1:p1);yy(1:p1);yy(p1+1:2*p1);yy(p1+1:2*p1)]; 
         ss((i-1)*bs+1:i*bs)=[is1;is1;is2];  
    end
    pse=sparse(ii,jj,ss,2*p1*m,2*p1*m);
    
end

end

 
