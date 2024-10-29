% Jacobi and Gauss-Seidel methods for solving Ax = b

A = [2 -1 0;-1 2 -1;0 -1 2];      % initial matrix and vector for
b = [7 2 5]';                   % Ax = b
k = 200;                       % number of iteration steps (if convergent)
x = [0 0 0]';            % initial guess (if convergent)
eps = 1*10^-8;           % to eliminate case e-val = 1
y = -1;            % For Jacobi, set y = 1. For Gauss-Seidel, set y = -1


n = size(A,1);
S = zeros(n,n);
if y == 1
    for i = 1:n
        S(i,i) = A(i,i);       % diagonal part of A
    end
    T = S - A;
elseif y == -1
    for i = 1:n
        for j = 1:i
            S(i,j) = A(i,j);    % lower triangular part of A
        end
    end
    T = S - A;
end
   
M = inv(S);
N = M*T;
v = eig(N);          % eigenvalues of inv(S)*T
w = -1;
for i = 1:size(v,1)
    if abs(v(i)) >= 1-eps       % check if mod of e-vals larger than 1
        if y == 1
            'Jacobi does not converge'
            w = 1;
            break
        elseif y == -1
            'Gauss-Seidel does not converge'
            w = 1;
            break
        end
    end
end

if w == -1            % if convergent, do iteration k times 
    for i = 1:k
        x = S\(T*x + b);  % solves linear system Sx = (T*x + b)  (better than finding inverse of S) (see matlab help)
    end
end

% end 

