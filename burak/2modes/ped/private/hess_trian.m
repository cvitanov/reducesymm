function [G,Q]=hess_trian(G)
% [G,Q]=hess_trian(G) transform a sequence of matrices
%into hessenberg-triangular form. 
%input: G
%       dimension of G: N x M*N. In the order: G1,G2,..Gm
%output: H, Q
%        both have the dimension: N x M*N. donated as 
%        H1,H2,..Hm. Q1,Q2,..Qm.
%       They satisfy:(Qi)^{T}*Gi*Q(i-1)=Hi. (Qm=Q0).
%       Therefore, (Qm)^{T}*G*Q0=H. where G=Gm*...G2*G1.

    N=size(G,1); M=size(G,2)/N; Q=repmat(eye(N),1,M); 
    for i=1:N-1
        for j=1:M-1 % householder transform on the ith column of G1,G2,..G(m-1). 
            [b,P]=househ(G(:,(j-1)*N+1:j*N),i,1); G(:,(j-1)*N+1:j*N)=b;
            G(:,j*N+1:(j+1)*N)=G(:,j*N+1:(j+1)*N)*P';
            Q(:,(j-1)*N+1:j*N)=Q(:,(j-1)*N+1:j*N)*P';  %update the orthogonal matrix.          
        end
        if i<N-1 %the second type householder transfrom on the ith column of Gm
            [b,P]=househ(G(:,(M-1)*N+1:M*N),i,2); G(:,(M-1)*N+1:M*N)=b;
            G(:,1:N)=G(:,1:N)*P'; Q(:,(M-1)*N+1:M*N)=Q(:,(M-1)*N+1:M*N)*P';
        end
    end 
    
    % deal with the roundoff errors to keep the Hess-trian structure.
    for j=1:M-1, G(:,(j-1)*N+1:j*N)=triu(G(:,(j-1)*N+1:j*N)); end
    G(:,(M-1)*N+1:M*N)=triu(G(:,(M-1)*N+1:M*N),-1);


end