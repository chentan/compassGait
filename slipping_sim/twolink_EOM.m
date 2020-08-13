function [D,C,G,B,H,dH,LfH,dLfH,J_com,J_hip] = twolink_EOM(x,params)
% Gives the matricies for the equations of motion
% of the two link walker
% Inputs:
% x - vector of positions and velocites, [q1,q2,q1',q2']
%
% Outputs
% D - mass / inertia matrix
% C - velocity
% G - Gravity
% B - torque transmission
% H - output functions
% LfH -
% dLfH -
%
% Added on 10 Aug 2009 for candidacy questions
% J_com - jacobian to give cartesian velocity of com
% J_hip - jacobian to give cartesian velocity of hip

%changed third input to params, structure of params, needed for disturbance
% testing
%
% dcp 13 may 2009
% recover original input
%
% 01 dec 2009
% Updated to also simulate the Simplest Walking model
% added several parameters to
% add in ground slope
% and manipulate equations correctly
% see SimplestWalker.nb
% will need beta1 and beta2
% because ignore the mass in one eqn and leave it in other
% See Garcia, et. al. 1998
%
% see twolink_EOM_old for previous working version
% also modified
% parameters()
%
% need to check impact map implications.

a = params.a;

[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);
U1 = x(3);
U2 = x(4);

qs = x(1:2);
qds= x(3:4);


% Commented rows are for changing from identical, symmetric legs
% changed so I could try simple walker in here
% parameters(); should work w/ the old way now

D = zeros(2);

D(1,1) = J + M*(L-LC)^2;
D(1,2) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
D(2,1) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
D(2,2) = 2*J + M*LC^2 + M*(L^2+(L-LC)^2-2*L*(L-LC)*cos(Q1));

% D(1,1) = J + m*(L-LC2)^2;
% D(1,2) = -J - (L-LC2)^2*m + L*(L-LC2)*m*cos(Q1);
% D(2,1) = -J - (L-LC2)^2*m + L*(L-LC2)*m*cos(Q1);
% D(2,2) = 2*J + 2*L^2*m - 2*L*LC2*m + LC2^2*m + LC^2*M - 2*L*(L-LC2)*m*cos(Q1);

G = zeros(2,1);
G(1) = GRAV*M*(L-LC)*sin(Q1-Q2);
G(2) = -GRAV*M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));

% G(1) = GRAV*m*(L-LC2)*sin(Q1-Q2-gam);
% G(2) = -GRAV*(L-LC2)*m*sin(Q1-Q2-gam) - GRAV*(L*m+LC*M)*sin(Q2+gam);%-GRAV*M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));

C = zeros(2);

C(1,2) = -L*M*(L-LC)*sin(Q1)*U2;
C(2,1) = -L*M*(L-LC)*sin(Q1)*(U1-U2);
C(2,2) = L*M*(L-LC)*sin(Q1)*U1;

% C(1,2) = L*(-L+LC2)*m*U2*sin(Q1); % -L*M*(L-LC)*sin(Q1)*U2;
% C(2,1) = -L*(L-LC2)*m*(U1 - U2)*sin(Q1) ;%-L*M*(L-LC)*sin(Q1)*(U1-U2);
% C(2,2) = L*(L-LC2)*m*U1*sin(Q1) ;%L*M*(L-LC)*sin(Q1)*U1;

B = [1; 0];


J_com = zeros(2,2);
J_com(1,1) = .5*(L-LC)*cos(Q1-Q2);
J_com(1,2) = .5*(L+LC)*cos(Q2) - .5*(L-LC)*cos(Q1-Q2);
J_com(2,1) = .5*(L-LC)*sin(Q1-Q2);
J_com(2,2) = -.5*(L+LC)*sin(Q2) - .5*(L-LC)*sin(Q1-Q2);
% J_com(2,1) = -.5*(L-LC)*sin(Q1-Q2);
% J_com(2,2) = -.5*(L+LC)*sin(Q2) + .5*(L-LC)*sin(Q1-Q2);

J_hip=zeros(2,2);
J_hip(1,2) = L*cos(Q2);
J_hip(2,2) = -L*sin(Q2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bezier Polynomials from Partial-Feedback Linearization

% Assign to variables for easier calculations
a0 = a(1);
a1 = a(2);
a2 = a(3);
a3 = a(4);
a4 = a(5);


% Convert from Bezier_Poly(x;a) to regular poly
% i.e. b0 + b1*x + b2*x^2 + ...
% Changed to correspond to a, dcp 21 APR 2009
% was hard-coded before
qp = -a4/2; % Angle of q2 at the beginning of a step
qm = -qp; % Angle of q2 at the end of a step
% Conversion to b0 + b1*x + b2*x^2 +b3*x^3 +b4*x^4
% Use the formula 6.7 Grizzle to expand the Bezier Poly
% It's given in a vector by 5 elements
b = [(qm+(-1).*qp).^(-4).*(a0.*qm.^4+qp.*((-4).*a1.*qm.^3+qp.*(6.*a2.* ...
    qm.^2+(-4).*a3.*qm.*qp+a4.*qp.^2))),(-4).*(qm+(-1).*qp).^(-4).*( ...
    a0.*qm.^3+(-1).*a1.*qm.^2.*(qm+3.*qp)+qp.*(3.*a2.*qm.*(qm+qp)+qp.* ...
    (a4.*qp+(-1).*a3.*(3.*qm+qp)))),6.*(qm+(-1).*qp).^(-4).*(a0.* ...
    qm.^2+a2.*qm.^2+4.*a2.*qm.*qp+(-2).*a3.*qm.*qp+a2.*qp.^2+(-2).* ...
    a3.*qp.^2+a4.*qp.^2+(-2).*a1.*qm.*(qm+qp)),(-4).*(qm+(-1).*qp).^( ...
    -4).*(a0.*qm+3.*a2.*qm+(-1).*a3.*qm+3.*a2.*qp+(-3).*a3.*qp+a4.*qp+ ...
    (-1).*a1.*(3.*qm+qp)),(a0+(-4).*a1+6.*a2+(-4).*a3+a4).*(qm+(-1).* ...
    qp).^(-4)];

% Again for easier calcs
b0 = b(1);
b1 = b(2);
b2 = b(3);
b3 = b(4);
b4 = b(5);
% Possibly make these lines into a subfunction or separate function

% Calculate outputs H to zero
H = x(1)+(-b0-b1*x(2)-b2*x(2)^2-b3*x(2)^3-b4*x(2)^4);

%look at H at contact

dH = zeros(1,4);
dH(1,2) = -(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3);
dH(1,1) = 1;

LfH = x(3)-(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3)*x(4);

dLfH = zeros(1,4);
dLfH(1,2) =  -(2*b2+6*b3*x(2)+12*b4*x(2)^2)*x(4);
dLfH(1,4) =  -(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3);
dLfH(1,3) = 1;
