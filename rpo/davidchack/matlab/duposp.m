function c = duposp(n,i);

%  npm = 1;  for i1 = 2:n, npm = npm*i1; end,
  npm = prod(1:n); 
  nsm = 2.^n;  per = 1:n;  sig = zeros(n,1);
  ii = mod(i,npm*nsm); ns = mod(ii,nsm); np = fix(ii/nsm); 

  for ip = n-1:-1:1,
    npm = npm/(ip+1);  kp = fix(np/npm);
    if kp > 0,   
      np = np - kp*npm; pt = per(ip-kp+1);
      for i1 = ip+(-kp+2:1),  per(i1-1) = per(i1);  end,
      per(ip+1) = pt;
    end,
  end,

  for is = n-1:-1:0,
    nsm = 2.^is;  ks = fix(ns/nsm); 
    if ks > 0,
      ns = ns - ks*nsm; sig(is+1) = -1;
    else
      sig(is+1) = 1;
    end,
  end,
  

  c = eye(n); c = c(per,:); c = diag(sig)*c;