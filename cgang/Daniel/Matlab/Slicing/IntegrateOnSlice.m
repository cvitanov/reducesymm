clear all
clc

figure(1)

% Set Porter-Knobloch system parameters and slice parameters
% ------------------------------------------------------------------
mu1 = -2.8023;
mu2 = 1; 
a1 = -1;
a2 = -2.6636;
b1 = 0; 
b2 = 0; 
c1 = -4.1122; 
c2 = 1.8854; 
e2 = 1;

% Set template 1
template1 = [1,1,0,0]'; 

%TODO: Set template 2
template2 = 3*randn(4,1);

% Set template_switching to 1 to enable two slice charting
template_switching = 0;

T = TPK; % Generator of infinitesimal rotations

% TODO: Rotate template 2 into slice of template 1
% [THETA,fval,template2] = PPShiftToSlice(0,template2,template1);
% template2
% 
% if (T*template2(:))'*(T*template1(:)) < 0 % Check directional condition on slice, if failed look for different plane crossing
%    [THETA,fval,template2] = PPShiftToSlice(pi,template2,template1);
% end

starting_template = 1;

params = {mu1, mu2, a1, a2, b1, b2, c1, c2, e2, template1, template2, starting_template};

disp('templates')
template1 
template2 

a = template2(:)'*(T*template1(:))
c = (T*template2(:))'*(T*template1(:))



% Set up initial condition and max integration time
% -------------------------------------------------------------------------
xinit = randn(4,1);

[THETA,fval,Xi] = PPShiftToSlice(0,xinit,template1); % Shift initial condition into slice of template 1



if (T*Xi)'*(T*template1) < 0 % Check directional condition on slice, if failed look for different plane crossing
    [THETA,fval,Xi] = PPShiftToSlice(pi,xinit,template1);
end

format long

% Set up integration    
X0 = [Xi; 0]' % Setup intial conditions including extra dimension for Phi
options = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events',@PKevents); % Setup integrator options

tstart = 0;
tfinal = 250;
tout = tstart;
xout = X0;
ieout = [];
teout = [];


 while tstart < tfinal  
     % Integrate one of the stop conditions (template change or crossing of slice border
     % or the maximum integration time is reached.
     [t,x,te,xe,ie] = ode45(@PKSliceEOM,[tstart tfinal],X0,options,params);
     
     % If stop condition is satisfied store integration up to that point
     % and check which stop condition has been reached
     nt = length(t); % Number of time steps taken between previous stop and current stop
     tout = [tout; t(2:nt)]; % Store time information up to stop event
     xout = [xout; x(2:nt,:)]; % Store trajectory
     ieout = [ieout; ie]; % Record which stop condition has been reached
     teout = [teout; te];
    %Display appropriate info for which stop condition was reached
    if ieout(end) == 1
        
        params{12} = 1;
        
        if tout(end) > 0 && tout(end) < (tfinal -.5)
            plot3(x(1:4:end,1),x(1:4:end,3),x(1:4:end,4),'r-'),hold on
            drawnow
        end
  
    elseif ieout(end) == 2
        
        if template_switching == 0
            params{12} = 1; 
        elseif template_switching == 1
            params{12} = 2;
        else
            error('template_switching is not a logical')
        end
        
        if tout(end) > 0 && tout(end) < (tfinal -.5)
            plot3(x(1:4:end,1),x(1:4:end,3),x(1:4:end,4),'b-'),hold on
            drawnow
        end
        
    elseif ieout(end) == 3
        disp('phidot blew up')
        length(tout)
    end % end if
    
    X0 = x(nt,:);    % Set conditions for when integration resumes
    tstart = t(nt);
    

 end % end while

% Plot Equilibria

for j = 1000:10:length(tout)
   phidot(j) = (xout(j,5)-xout((j-1),5))/(tout(j)-tout(j-1));
%     if abs(phidot(j)) > 5
%         ,subplot(1,2,1)
%         plot3(xout(j,1),xout(j,2),xout(j,5),'rx')
%         
%         ,subplot(1,2,2)
%         plot3(xout(j,3),xout(j,4),xout(j,5),'rx')
%     end
%     hold on
end
    
   
