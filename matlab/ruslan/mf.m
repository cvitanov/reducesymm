function [aa, th] =mf(a, slc)
%
% slc: Fourier mode to use for fixing a slice by c_slc=0.
ac=a(1:2:size(a))+1i*a(2:2:size(a));
phi=angle(ac(slc))/slc;
%theta=-atan2(a(2*slc),a(2*(slc-1)+1)); %theta=-phi
as = exp(-1i*(1:size(ac))'*phi).*ac;
aa=zeros(size(a));
aa(1:2:size(a))=real(as);
aa(2:2:size(a))=imag(as);
if nargout > 1,
    th=phi;
end
end
