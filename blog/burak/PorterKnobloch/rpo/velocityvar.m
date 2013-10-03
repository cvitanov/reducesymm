function v = velocityvar(x)
%Variational velocity function
%Relevant notes: %http://www.cs.colorado.edu/~lizb/chaos/variational-notes.pdf

n=4; % System dimension, change accordingly

v = velocity(x(1:n));
A = StabilityMatrix(x(1:n));

VariationsMatrix = vector2matrix(x(n+1:size(x,1)));
MatrixVelocity = A*VariationsMatrix;
n2Velocity = matrix2vector(MatrixVelocity);

v(n+1:n^2+n) = n2Velocity;
