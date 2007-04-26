  function tall(f);
%	TALL(F) - make the axis taller by a fraction F

  ax = get(gca,'pos');
  ax(2) = ax(2) + (1-f)*ax(4)/2;
  ax(4) = ax(4)*f;
  set(gca,'pos',ax);