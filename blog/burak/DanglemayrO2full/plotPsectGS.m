function void = plotPsectGS(filename, V);
%Projects the Poincare Section coordinates onto the basis stored in V,
%and then plots the first two. Hence V must have the perpendicular 
%coordinate in the third row.

warning ("off", "Octave:broadcast"); % This function uses broadcasting purposefully
load(filename); % Load the file with the data set

psect = ps(1:3,:); % Seperate time and position entries
psectrelative = psect(1:3,:)-xhatp; %Positions measured relative to the temp.

psectprojected(1:3,:) = V' * psectrelative(1:3,:);

plot(psectprojected(1,:),psectprojected(2,:), '.', 'markersize', 5)

void = 1;
