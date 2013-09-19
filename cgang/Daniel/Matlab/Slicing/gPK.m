% This generates the rotation matrix for a finite angle theta
% about the z axis for the Porter Knobloch system in the 
% (x1,x2,y1,y2,z) basis. Works for both symbolic and numeric
% inputs.

function g = gPK(theta)

    g(1,1) = cos(theta);
    g(1,2) = sin(theta);
    g(1,3) = 0;
    g(1,4) = 0;
    
    g(2,1) = -sin(theta);
    g(2,2) = cos(theta);
    g(2,3) = 0;
    g(2,4) = 0;
   
    g(3,1) = 0;
    g(3,2) = 0;
    g(3,3) = cos(2*theta);
    g(3,4) = sin(2*theta);
    
    g(4,1) = 0;
    g(4,2) = 0;
    g(4,3) = -sin(2*theta);
    g(4,4) = cos(2*theta);
end