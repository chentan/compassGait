% MAIN ACROBOT RUNNING SCRIPT (run this to execute simulation)
% A GENERAL COMPASS GAIT BIPED THAT ALLOWS FOOTSLIPPING
%
% These parameters work for point-foot. HZD controller is adopted. For more
% details about HZD controller, refer to "Feedback control of dynamic
% bipedal robot locomotion"
%
% x(1) = swing leg, x(2) = stance leg, x(3) = x1dot, x(4) = x2dot
% Desired step length = 0.4450

close all
clear
clc

%% set up friction-related parameters
% scale_coeff: ration between static and kinetic friction coefficients
params.scale_coeff = 1.2;

% fric_coeff: coefficient of kinetic friction 
% params.fric_coeff = .35/1.2; % change ground friction HERE
params.fric_coeff = .4/1.2; % change ground friction HERE

params.vary_fric = 0; 
params.flags.foot_slip = 0; % start with no-slip

%%
% gait parameter
params.a_0 = [0 0 1.12 0.48 0.4488]; % this is a feasible gait on rough ground

% disturbance
params.vel.dist = 0; % disturbance magnitude
params.vel.dir = 1;  % what direction is disturbance aligned with (see v_disturbance)
params.vel.cycle_percent = 0.6; % where I want the disturbance to happen
params.flags.vel_dist = 0;    % existence of velocity disturbance? 1 Yes 0 No

params.flags.step_one = 1 ; 
 
% Run the walker
step_num = 5; % change number of steps here
data = run_walker4(params.a_0,step_num,params); 
s1 = "Success";
s2 = 'Failure-FallBack';
s3 = 'NegativeFn';

if strcmp(data.comment,s1)
    disp('Success');
elseif strcmp(data.comment,s2)
    disp('Fall backward');
elseif strcmp(data.comment,s3)
    disp('Negative force required')
else
    disp('N/A')
end
 
% Animation
animate3(data);