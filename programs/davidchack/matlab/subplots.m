function hax = subplots(nr, nc, be, bi)
% SUBPLOTS Create axes in tiled positions with specified external and
%          internal boundaries.
%   HAX = SUBPLOTS(NR, NC, BE, BI) creates NR-by-NC matrix of axes with
%   BE = [BEL BER BEB BET] (left, right, bottom, top) boundaries around
%   all axes and BI = [BIL BIR BIB BIT] around each axis.
%   HAX is a vector of handles to the axes.

% Ruslan L. Davidchack,   Apr. 28, 2001 

  wid = (1.0 - be(1) - be(2))./nc - bi(1) - bi(2);
  hei = (1.0 - be(3) - be(4))./nr - bi(3) - bi(4);
  hax = zeros(nr.*nc,1);  cnt = 0;
  clf; 
  for ir = 1:nr,
    for ic = 1:nc,
      cnt = cnt + 1;
      left = be(1) + bi(1) + (ic-1).*(wid+bi(1)+bi(2));
      bot = be(3) + bi(3) + (nr-ir).*(hei+bi(3)+bi(4));
      hax(cnt) = axes('position',[left bot wid hei]);
  end, end,

return;
