function M = vector2matrix(v)

%This function converts an n^2 dim vector to an n x n matrix
%Conversion is consistent with the notes given here:
%http://www.cs.colorado.edu/~lizb/chaos/variational-notes.pdf

n=length(v);
n=sqrt(n);

for i = 1:n
	for j = 1:n

	M(i,j) = v((j-1)*n+i);
	
	end
end
