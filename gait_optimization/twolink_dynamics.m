function dx = twolink_dynamics(t,x,params)
% dx = twolink_dynamics(t,x,params)
% dx = numerical derivatives
% t = time (so far unused in calculations)
% x = states
% params = structure of model and control parameters 
% (see ...)

% dcp 03 FEB 2010 
% comment code so Anne can use to simulate her impact maps and etc.
% 

% changed third input to params, structure of params, needed for disturbance
% testing
%
% dcp 13 may 2009
% recover original input
% a = params.a;

%commented for running optimizations
% global t_2 torque 
% global y forces

%% Load EOM entries needed to compute
% dx = f(x) + g(x)*u
params.a;
[D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(x,params);
% [D,C,G,B,H,~,LfH,dLfH] = twolink_EOM(x,params);
% 2nd line only compatible in newer versions of matlab

Fx = D\(-C*x(3:4)-G); % faster than inv(D)*
% Fx = inv(D)*(-C*x(3:4)-G);

Gx = D\B; % faster than inv(D)*
% Gx = inv(D)*B;

% Calculate Lie Derivatives needed for control
LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];
%% Calculate control input
% Alternate control ideas are commented

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Using PD Controller
% % % try eqn 5.97
% % Ernie uses low-level PD control 
 Kd = 25;
%  Kd = 50; % Tan changed
 Kp = 25; % JN, changing this will decrease tracking errors but drive up control inputs
% Kp = 250; % Tan changed
 e=.2;
% e=.02; % Tan changed
 v = pd_control(Kp,Kd,e,H,LfH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bernstein-Bhat controller (uses feedback linearization)
% Copied / Modded from Eric's 3-link
% This performs better than PD in simulation
% Not sure if its actually used on the hardware
% Use this to match answers in Eric's book
%[v]=control_two_link(H,LfH);
u = inv(LgLfH)*(v-Lf2H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Use a switching controller
% % turns on or off the CA controller, based on criterion
% % with no flag, should run regular Bernstein-Bhat control
% %
% u = switching_control(x,params);
% % Need this check in case my switching rules suck
% if isnan(u)
%     disp('woah!')
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Hack for passive walker
% u = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute dx
dx(1:2) = x(3:4);
dx(3:4) = Fx+Gx*u;
%% Mod to calc angular momentum and store as state 5

[L,LC,M,J,GRAV] = parameters();
% [L,LC,M,~,~] = parameters(); %Again can be used in newer matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Use these if trying to run simple walker
% [L,LC1,LC2,M,m,J,GRAV,gam] = parameters(); 
% LC = LC1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q1=x(1);
Q2=x(2);
Q1d = x(3);
Q2d = x(4);

dtemp = [(1/2).*(L.^2+(-1).*LC.^2).*M.*(Q1d+(-2).*Q2d).*sin(Q1),0,(-1/2).* ...
  M.*((-1).*L.^2+2.*L.*LC+(-1).*LC.^2+(L.^2+(-1).*LC.^2).*cos(Q1)),( ...
  -1/2).*M.*(2.*L.^2+2.*LC.^2+(-2).*(L.^2+(-1).*LC.^2).*cos(Q1))];

% Rate of change of angular momentum (dot ang mom)
dangM = dtemp*dx';
% Alternate way to calc
% dangM = dtemp*[x(3:4);Fx] + dtemp*[zeros(2,1);Gx]*u; 
dx(5) = dangM;

% Eventually make a function of of this, such that the ang mom is based on
% the model and not hard coded here
%% Reorient
dx=dx.';

%% Old way of tracking control signal and force
% This code is probably not needed now
% but make sure before deleting
% uncomment later 
% torque = [torque ; u.'];
% t_2 = [t_2 ; t];
% y = [y ; H.'];
% [FT,FN] = stance_force_two_link(x,dx);
% % FT=x(5);
% % FN=x(6);
% forces = [forces ; FT FN];