  function wide(f);
%	WIDE(F) - make the axis wider by a fraction F

  ax = get(gca,'pos');
  ax(1) = ax(1) + (1-f)*ax(3)/2;
  ax(3) = ax(3)*f;
  set(gca,'pos',ax);