function void = plotPsectGS(filename, V);
%Projects the Poincare Section coordinates onto the basis stored in V,
%and then plots the first two. Hence V must have the perpendicular 
%coordinate in the third row.

load(filename); % Load the file with the data set

psect = ps(1:3,:);

psectprojected(1:3,:) = V' * psect(1:3,:);

plot(psectprojected(1,:),psectprojected(2,:), '.', 'markersize', 5)
title('Poincare section')

void = 1;
