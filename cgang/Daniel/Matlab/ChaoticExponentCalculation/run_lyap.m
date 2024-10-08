% run_lyap.m
% by Daniel Borrero
% ------------------------------------------------------------
% This script calculates the Lyapunov exponents for the 4-D
% Porter-Knobloch system and helps optimize parameters for chaos by
% varying them and checking which direction the Lyapunov exponents move.
% Below is a mini-blog of how I got to a set of chaotic parameters starting
% from a candidate set that was listed in the Chaos gang blog.

% Starting point 8/21/2013 10:36 AM - Trying to make 3rd and 4th lyapunov
% exponents more negative while keeping the leading exponent positive and
% the second exponent close to zero (within 0.002). Noise amplitude for
% parameter variation is 0.1. Starting all calculations from [1 1 1 1]
% params = [-0.38 0.38 -1.31 -2.60 1.504 0.22 -1.00 1.00 1.6100];
% lyapunov_exp =[0.0428    0.0011   -0.0091   -0.2789];

% Achieved parameter set that gives comparable Lyapunov exponents to the
% Rossler system (0.09,0,-9.77).

% Starting point 8/21/2013 7:45 PM - Trying to increase the magnitude of the
% leading exponent while keeping the second exponent within 0.001 of zero
% and the third exponent more negative than -8. 
% params = [-2.8023    0.6128   -0.9146   -2.6636   -0.6141   -0.0144   -4.1122    1.8854    2.1769];
% lyapunov_exp = [0.0928    0.0004   -9.1393  -17.7721];

% Didn't go anywhere after 1400 iterations... 

% Starting point 8/22/2013 7:39 AM - Same as before. Increasing noise 
% amplitude to 0.3 to see if I can get out of this local minimum.

% Hmmmm... still didn't change after 600 iterations...

% Starting point 8/22/2013 10:28 AM - Same as before. Increasing noise 
% amplitude to 0.5 to see if I can get out of this local minimum.

% No luck, but already getting pretty chaotic looking dynamics.

% params = [mu1 mu2 a1 a2 b1 b2 c1 c2 e2];
load params.mat
load lyapunov_exp.mat

warning off

params(10) = 1; 

for i = 1:100
    testparams = params ;%+ 0.5*(2*rand(1,9)-1);
    u = 3*rand(1);
    v = 3*rand(1);
    phi = pi/2*rand(1);
    w = 2*u*sqrt(v)*cos(phi);
    q = 2*u*sqrt(v)*sin(phi);
    [T,Res]=lyapunov(4,@PKeqInv,@ode45,0,0.5,2500,[u v w q],0, testparams);
   
    Res(end,:)
    
    L = [L Res(end,1)]
    
%     if  Res(end,1) > lyapunov_exp(1) && Res(end,3) < -8 && abs(Res(end,2)) < 0.001
%         disp('Win')
%         params = testparams
%         lyapunov_exp = Res(end,:)
%         save params.mat params
%         save lyapunov_exp.mat lyapunov_exp
%         
%     end
end


