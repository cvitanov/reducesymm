function aa =fmRot(a, phi)

ac=a(1:2:size(a))+1i*a(2:2:size(a));
as = exp(1i*(1:size(ac))'*phi).*ac;
aa=zeros(size(a));
aa(1:2:size(a))=real(as);
aa(2:2:size(a))=imag(as);
end