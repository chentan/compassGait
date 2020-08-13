function dx = twolink_dynamics(t,x,params)
% dx = numerical derivatives
% t = time (so far unused in calculations)
% x = states
% dx = f(x) + g(x)*u

params.a;
[D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(x,params);

Fx = D\(-C*x(3:4)-G); 

Gx = D\B; 

% Calculate Lie Derivatives needed for control
LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];

%% Calculate control input
Kd = 25;
Kp = 25; 
e = .2;
v = pd_control(Kp,Kd,e,H,LfH);

% Bernstein-Bhat controller (uses feedback linearization)
% Copied / Modded from Eric's 3-link
% This performs better than PD in simulation
% Not sure if its actually used on the hardware
% Use this to match answers in Eric's book
%[v]=control_two_link(H,LfH);
u = inv(LgLfH)*(v-Lf2H);

%% Compute dx
dx(1:2) = x(3:4);
dx(3:4) = Fx+Gx*u;

%% Mod to calc angular momentum and store as state 5
% [L,LC,M,J,GRAV] = parameters();
% 
% Q1 = x(1);
% Q2 = x(2);
% Q1d = x(3);
% Q2d = x(4);
% 
% dtemp = [(1/2).*(L.^2+(-1).*LC.^2).*M.*(Q1d+(-2).*Q2d).*sin(Q1),0,(-1/2).* ...
%     M.*((-1).*L.^2+2.*L.*LC+(-1).*LC.^2+(L.^2+(-1).*LC.^2).*cos(Q1)),( ...
%     -1/2).*M.*(2.*L.^2+2.*LC.^2+(-2).*(L.^2+(-1).*LC.^2).*cos(Q1))];

% % Rate of change of angular momentum (dot ang mom)
% dangM = dtemp*dx';
% dx(5) = dangM;

%% Reorient
dx=dx';