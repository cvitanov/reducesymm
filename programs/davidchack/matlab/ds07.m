% For DS07 talk
%% Troughness of |F| for Ikeda map
 clear;  par = [1.0; 0.9; 0.4; 6.0];
 bndr = [-1 1.5 -2 1];
 [x0,y0] = bndrgrid(bndr,160000);
 for p = 4:4:12,
   [xx,yy] = ikedaxy(par,p,x0,y0);
   figure(1); clf;
   F = (xx-x0).^2+(yy-y0).^2;
   pcolor(x0,y0,F); caxis([0 5]); 
   shading flat; axis(bndr); pixgrid(size(x0),0);
   print(['ik16p' num2str(p) 'norm.emf'],'-dmeta'); pause;
 end,
