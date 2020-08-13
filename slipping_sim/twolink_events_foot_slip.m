function [value,isterminal,direction] = twolink_events_foot_slip(t,x,params)
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
    if abs(x(1)) < 0.05
        isterminal(1) = 0; % don't stop
    else
        isterminal(1) = 1; % can stop integration
    end
end
direction = -1;

%% hips get too close to ground--kill simulation
value(2) = L*cos(x(2))-0.5*L;
isterminal(2) = 1;
direction(2) = -1;

%% kill because falling backwards
value(3) = x(2) + a(5)/2 + 0.01; 
isterminal(3) = 1;
direction(3) = -1;

%% distrubance
TH_val = -a(5)/2+params.vel.cycle_percent*a(5);
if params.flags.vel_dist == 1
    value(4) = x(2)-TH_val;
    isterminal(4) = 1;
    direction(4) = 1;
else
    value(4) = 1;
    isterminal(4) = 0;
    direction(4) = 1;
end

%% Check whether the slipping foot velocity is zero
% Stop slip only when the velocity is zero and the required constraint
% force can be met by the friction
% To tackle the numerical error, assume velocity at 1e-9 is just zero.
% value(5) = x(6)-sign(x(6))*1e-9;
value(5) = abs(x(6))-1e-9;
x_check = [x(1:2)',x(4:5)'];
dx_check = twolink_dynamics(0,x_check',params);
[FT_init,FN_init] = stance_force_two_link(x_check,dx_check);
% % Note that \miu_s > \miu_k
% if abs(FT_init/FN_init) > (coeff_of_fric)*scale_coeff
%     isterminal(5) = 0; % don't stop integration
% else
%     isterminal(5) = 1; % can stop integration
% end
% switch to No Slip when detecting slipping velocity is zero
isterminal(5) = 1;
% direction(5) = 0;
direction(5) = -1;

%% Integration issue: too large joint velocity
% if x(4) ~= 0
%     value(6)= x(4)-sign(x(4))*1e2;
%     isterminal(6) = 1;
%     direction(6) = 0;
% else
%     isterminal(6) = 0;
%     direction(6) = 0;
% end
value(6)= abs(x(4))-1e2;
isterminal(6) = 1;
direction(6) = 1;