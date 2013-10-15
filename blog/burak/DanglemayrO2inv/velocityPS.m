function v = velocityPS(x,i)
%Transformed velocity function according to ChaosBook eq. 3.10
%ith position variable is choosen as the time variable
%State vector: xnew = [x1, x2, ..., t, ... ,xd]
%dt/dxi = 1/vi
%dxd/dxi = vd / vi

v = velocity(x);
vi = v(i);
v = v./vi;
v(i) = 1/vi;
