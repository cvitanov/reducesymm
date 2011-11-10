function [ts, as] = shiftrpo(tt, aa, L, tend, ph, ishift)
%  Shift initial point of RPO to tt(ishift).
  if ishift > 1 & ishift <= length(tt),
    ts = [tt(ishift:end) tt(2:ishift)+tend]-tt(ishift);
    aas = aa(:,2:ishift);
    vs = (aas(1:2:end,:)+1i*aas(2:2:end,:));
    vs = vs.*repmat(exp((2i*pi/L).*ph.*(1:size(vs,1))'),1,ishift-1);
    aas(1:2:end,:) = real(vs);  aas(2:2:end,:) = imag(vs);
    as = [aa(:,ishift:end) aas]; end