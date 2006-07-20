function C = nearsd(Q);
%  Find Schmelcher-Diakonos (SD) matrix nearest to the orthogonal matrix Q

% (c) Ruslan L Davidchack, 02-Mar-2005

  C = Q;  n = size(Q,1);  sp = zeros(n,2);
  for ii = 1:n,
    [m,im] = max(abs(C(:)));
    ic = floor((im-1)/n)+1;  ir = mod((im-1),n)+1;
    sp(ir,1) = ic;  sp(ir,2) = sign(C(ir,ic));
    C(ir,:) = 0;  C(:,ic) = 0;
  end,
  C((sp(:,1)-1)*n+(1:n)') = sp(:,2);