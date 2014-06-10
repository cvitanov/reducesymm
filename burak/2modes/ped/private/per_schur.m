function [M,Q]=per_schur(M,Q,maxnp,eps)
% function [M,Q]=xdper_schur(MM,QQ,L,U) transforms an unreduced
% hessenberg-triangular sequence of matrices into periodic Schur form.
% This iteration method is based on the Implicit-Q theorem.
%input: 
%       MM: unreduced hessenberg-triangular sequence MM1,MM2,MM3..MMm. the
%           first m-1 matrices are upper-triangular, and the last one is
%           hessenberg matrix.
%       QQ: the initial orthogonal matrices.
%       maxnp: maximal number of iterations, default value is 100.
%       eps: the tolerance small number.
%output:
%       M: the periodic Schur form. M1,M2,..M(m-1) are upper-triangular
%       matrices.
%       Q: the orthogonal matrices, such that,
%       Mi=Qi^{T}*QQi*MMi*QQ(i-1)^{T}*Q(i-1)
 n=size(M,1);
 if nargin<4, eps=10^-20; end
 if nargin<3, maxnp=1000; end
[M,Q,~]=xdper_schur_onestep(M,Q,1,n,maxnp,eps);
end
%---------------------------------------------------------------------------
function [M,Q,nit]=xdper_schur_onestep(M,Q,L,U,maxnp,eps)
%       L : the start index.
%       U : the end index.
%Initialize.
n=size(M,1); m=size(M,2)/n; nit=0;
%deal with the dimension==1 case.
if U-L==0,  return, end
%deal with the case: dimension=2. if the eigenvalues are complex pairs,
%no further reduction is needed; otherwise, we need to turn it into 
%diagonal form.
if U-L==1
    mhess=eye(2);
    %normalize the matrix to avoid overflow or downflow
    for j=1:m, mhess=M(L:U,(j-1)*n+L:(j-1)*n+U)*mhess; mhess=mhess/max(max(abs(mhess))); end
    del=(mhess(1,1)-mhess(2,2))^2+4*mhess(1,2)*mhess(2,1); 
    %judge the eigenvalues are complex pair or real numbers.
    if del <0, return, 
    else
        [ve,la]=eig(mhess);
        if la(1,1)> la(2,2), [R,~]=givens(ve(:,1),n,L,U); % let the larger eigenvelue be in the upper place.
        else [R,~]=givens(ve(:,2),n,L,U);
        end
        M(:,(m-1)*n+1:m*n)=R*M(:,(m-1)*n+1:m*n);
        M(:,1:n)=M(:,1:n)*R'; Q(:,(m-1)*n+1:m*n)=Q(:,(m-1)*n+1:m*n)*R';
        for k=1:m-1 % restore the hessenberg-triangular form
            [R,y]=givens(M(L:U,(k-1)*n+L),n,L,U); M(:,(k-1)*n+1:k*n)=R*M(:,(k-1)*n+1:k*n);
            M(L:U,(k-1)*n+L)=y; % elimilate roundoff errors.
            M(:,k*n+1:(k+1)*n)=M(:,k*n+1:(k+1)*n)*R'; Q(:,(k-1)*n+1:k*n)=Q(:,(k-1)*n+1:k*n)*R';
        end
        %Afther the similarity transform. M(U,(m-1)*n+L) should be
        %zero,but round error exits => cannot return now. 
        M=formhm(M,eps); 
        %fprintf(1,'2x2 has real eigenvalues: %d\t%d\n', L,U);       
    end    
end 

id=L-1; np=0;
while size(id,2)<3 && np<maxnp
    %Here define the shift. Right now don't use any shift.    
    H=M(L:L+1,(m-1)*n+L);    
    [R,y]=givens(H,n,L,L+1); M(:,(m-1)*n+1:m*n)=R*M(:,(m-1)*n+1:m*n);  M(L:L+1,(m-1)*n+L)=y;
    M(:,1:n)=M(:,1:n)*R'; Q(:,(m-1)*n+1:m*n)=Q(:,(m-1)*n+1:m*n)*R';    
   
    for j=L:U-1      
        for k=1:m-1
            [R,y]=givens(M(j:j+1,(k-1)*n+j),n,j,j+1); M(:,(k-1)*n+1:k*n)=R*M(:,(k-1)*n+1:k*n); M(j:j+1,(k-1)*n+j)=y;            
            M(:,k*n+1:(k+1)*n)=M(:,k*n+1:(k+1)*n)*R'; Q(:,(k-1)*n+1:k*n)=Q(:,(k-1)*n+1:k*n)*R'; 
        end
       
        if j<U-1
            % check whether M is already the hessenberg-triangular form. 
             if abs(M(j+2,(m-1)*n+j))< eps*(abs(M(j+1,(m-1)*n+j))+abs(M(j+2,(m-1)*n+j+1)))
                 M(j+2,(m-1)*n+j)=0; break, 
             end
            [R,y]=givens(M(j+1:j+2,(m-1)*n+j),n,j+1,j+2); M(:,((m-1)*n+1:m*n))=R*M(:,((m-1)*n+1:m*n)); M(j+1:j+2,(m-1)*n+j)=y;
            M(:,1:n)=M(:,1:n)*R'; Q(:,(m-1)*n+1:m*n)=Q(:,(m-1)*n+1:m*n)*R';   
        end  
     end
    %  get rid of round error.
    M=formhm(M,eps);
    
    id=L-1;
    for j=L:U-1, if M(j+1,(m-1)*n+j)==0, id=[id,j]; end, end, id=[id,U];
    %if size(id,2)>2, disp([L,U,size(id,2), id]);  end   
    if size(id,2)>2
        for j=1:size(id,2)-1, [M,Q,nitt]=xdper_schur_onestep(M,Q,id(j)+1,id(j+1),maxnp,eps); nit=nit+nitt; end
    end
    np=np+1; nit=nit+1; 
end
%fprintf(1,'nit=%d\t np=%d\t L=%d\t U=%d\n',nit,np,L,U);
end

%-------------------------------------------------------------------------
function M=formhm(M,eps)
% deal with the rundoff errors to keep the Hess-trian structure.
n=size(M,1); m=size(M,2)/n; if nargin<2,eps=10^-20; end
for j=1:m-1, M(:,(j-1)*n+1:j*n)=triu(M(:,(j-1)*n+1:j*n)); end
M(:,(m-1)*n+1:m*n)=triu(M(:,(m-1)*n+1:m*n),-1);
% the following will recover the zeros in the subdiagonal, which could be
% ruined by the round error introduce by function recursion.
for j=1:n-1
    if abs(M(j+1,(m-1)*n+j))< eps*(abs(M(j,(m-1)*n+j))+abs(M(j+1,(m-1)*n+j+1)))
       M(j+1,(m-1)*n+j)=0;          
    end
end
end
%--------------------------------------------------------------------------
function N=sform(M,n,flag)
    if flag==1
        N=triu(M);
        N(n+1,n)=M(n+1,n);
    else
        N=triu(M,-1);
        if n+2<=size(M,1), N(n+2,n)=M(n+2,n); end 
    end            
end

