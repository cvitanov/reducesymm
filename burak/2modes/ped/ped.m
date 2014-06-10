function [pereig,ve]=ped(J)

    flag=0;
    if nargout==2,flag=1;end
    
    J=reverse(J); % reverse the order of J sequence.
    [M,Q]=hess_trian(J); % transformation to Hessenberg-triangular form.
    [tM,tQ]=per_schur(M,Q,1000,10^-14); % periodic Schur form.

    if(flag==0) pereig=per_eigvalue(tM);
    else [pereig,ve]=per_eigreorder(tM,tQ); %reordering algorithm for eigenvectors.
    end

end
