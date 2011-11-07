function [aa] = ksfmRefl(a)
% Apply reflection operator to KS point in Fourier modes representation.

aa=a;
aa(1:2:end-1)=-a(1:2:end-1);

end

