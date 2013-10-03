function v = matrix2vector(M)

%This function converts an n x n matrix to an n^2 dim vector 
%Conversion is consistent with the notes given here:
%http://www.cs.colorado.edu/~lizb/chaos/variational-notes.pdf

n=size(M,1);

for i = 1:n
	for j = 1:n
	
	v((j-1)*n+i)=M(i,j);
	
	end
end

v=v';
