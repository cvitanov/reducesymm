function [R,y]=givens(x,N,i,j)
%[R,y]=givens(x,N,i,j): Givens rotation.
%Input: x: a 2*1 vector.
%       N: dimension of the matrix.    
%       i,j: the index of x in the original matrix.
%Output: R: the rotation matrix.
%       y: the rotated 2*1 vector.

% the consideration below is for overflow or downflow. so  verytime get the
% rotation matirx, the rotation cannot make sure that x(2)=0. you should 
% substitute x=y.
    %eps=10^-160;
    eps=0;
    if i>=j, disp('givens rotation err: i>=j'); end    
    if abs(x(2))<eps*abs(x(1)), R=eye(N); y=[x(1);0];
    elseif abs(x(1))<eps*abs(x(2)),
        R=eye(N);R(i,i)=0;R(j,j)=0; R(i,j)=sign(x(2));R(j,i)=-sign(x(2));
        y=[0;abs(x(2))];
    else
        R=eye(N);
        nor=sqrt(x(1)*x(1)+x(2)*x(2));
        R(i,i)=x(1)/nor; R(j,j)=x(1)/nor; R(i,j)=x(2)/nor; R(j,i)=-x(2)/nor;
        y=[nor; 0];
    end
end