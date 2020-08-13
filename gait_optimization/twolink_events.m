function [value,isterminal,direction] = twolink_events(t,x,params)
%changed third input to params, structure of params, needed for disturbance
% testing
%
% dcp 13 may 2009
% recover original input
a = params.a;

%update this to include new switching control parameters
[th3d,th2d,alpha,epsilon] = control_params_walker;
% 	[L,LC1,LC2,M,m,J,GRAV,gam] = parameters();
[L,LC,M,J,GRAV] = parameters();
[pFoot1,pFoot2,pHip,pTorso]=limb_pos(x,0);

%     value(1) = th2d - x(2);

% Check if swing foot has hit ground
value(1) = pFoot2(2);
% value(1) = double(abs(pFoot2(2))>1e-2); % Tan changed, too large acuracy
% value(1) = double(abs(pFoot2(2))>5e-6); % Tan changed, why not work?

% to prevent heel scuffing the simulation is only stop if the stance leg
% is past vertical
if x(2) < -0.001
    % if x(2) < 0.2 % Tan changed
    isterminal(1) = 0; % don't stop
else
    % I want to check whther abs(x(1))==0 here actually, but not convenient
    % to express it directly, since x(1)can never be perfectly zero. And
    % there is integration issue here.
    if abs(x(1)) < 0.05
        isterminal(1) = 0; % don't stop
    else
        isterminal(1) = 1; % can stop integration
    end
end

% hips get too close to ground--kill simulation
value(2) = L*cos(x(2))-0.5*L;
isterminal(2) = 1;
direction(2) = -1;


value(3) = x(2) + a(5)/2;  % kill because falling backwards
isterminal(3) = 1;
direction(3) = -1;

value(4:5)= [1 1];  % fix disappearing events problem?
isterminal(4:5) = [0 0];
direction(4:5) = [1 1];

%% Integration issue: too large joint velocity
% if x(3) ~= 0
%     value(6)= x(3)-sign(x(3))*1e2;
%     isterminal(6) = 1;
%     direction(6) = 0;
% else
%     value(6)= 1;
%     isterminal(6) = 0;
%     direction(6) = 0;
% end
value(6)= abs(x(3))-1e2;
isterminal(6) = 1;
direction(6) = 1;
end