% This function generates the generator of infinitesimal rotations
% for SO(2)for the Porter Knobloch system in the (x1,x2,y1,y2,z) basis
function T = TPK
T = zeros(4,4);

T(1,2) = 1;
T(2,1) = -1;
T(3,4) = 2;
T(4,3) = -2;

end

