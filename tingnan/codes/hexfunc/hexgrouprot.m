function [ vec_out] = hexgrouprot(grot, basevec)
% rotate and translate the vector according to the group action
% basevec should be a column vector.
% the cumulated grot will determine which translation action to take:
% grot = n , rot(basevec, pi/3*n);

theta = pi/3*grot;

rotmat = [cos(theta), -sin(theta);
          sin(theta), cos(theta)];

vec_out = rotmat * basevec;

end

