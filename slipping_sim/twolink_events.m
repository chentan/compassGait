function [value,isterminal,direction] = twolink_events(t,x,params)
coeff_of_fric = params.fric_coeff;
scale_coeff = params.scale_coeff;

a = params.a;

% update this to include new switching control parameters
[th3d,th2d,alpha,epsilon] = control_params_walker;
[L,LC,M,J,GRAV] = parameters();
[pFoot1,pFoot2,pHip,pTorso] = limb_pos(x,0);

%% Check if swing foot has hit ground
value(1) = pFoot2(2);
if x(2) < -0.001
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
% Stop when the foot hit the ground from upwards
direction = -1;

%% hips get too close to ground--kill simulation
value(2) = L*cos(x(2))-0.5*L;
isterminal(2) = 1;
direction(2) = -1;

%% kill because falling backwards
% +0.01 to make this event incompatible with event 1
value(3) = x(2) + a(5)/2 + 0.01;  
isterminal(3) = 1;
direction(3) = -1;

%% disturbance
TH_val = -a(5)/2 + params.vel.cycle_percent*a(5);
% added because problems w/ multiple disturbance per step
if params.flags.vel_dist==1
    % Velocity Disturbance
    % stop and adjust when pass the disturbance angle
    value(4) = x(2) - TH_val;
    isterminal(4) = 1;
    direction(4) = 1;
else
    value(4) = 1;
    isterminal(4) = 0;
    direction(4) = 1;
end

%% Check the required friction coefficient at each instant
% Works for the case from no-slip to slip
dxout = twolink_dynamics(0,x,params);
[FT,FN] = stance_force_two_link(x,dxout);
fric_coef = abs(FT/FN);
% Note that \miu_s > \miu_k
value(5) = fric_coef - (coeff_of_fric)*scale_coeff;
isterminal(5) = 1;
direction(5) = 1;

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