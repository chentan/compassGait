function [D,C,G,B1,B2,H,dH,LfH,dLfH,J_com,J_hip] = twolink_EOM_foot_slip(x,params)

a = params.a;

[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);
Q3 = x(3);
U1 = x(4);
U2 = x(5);
U3 = x(6);

qs = x(1:3);
qds= x(4:6);


% Commented rows are for changing from identical, symmetric legs
% changed so I could try simple walker in here
% parameters(); should work w/ the old way now

D = zeros(3);

D(1,1) = J + M*(L-LC)^2;
D(1,2) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
D(1,3) = M*(L-LC)*cos(Q1-Q2);
D(2,1) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
D(2,2) = 2*J + M*LC^2 + M*(L^2+(L-LC)^2-2*L*(L-LC)*cos(Q1));
D(2,3) = -M*(L-LC)*cos(Q1-Q2) + M*(L+LC)*cos(Q2);
D(3,1) = M*(L-LC)*cos(Q1-Q2);
D(3,2) = -M*(L-LC)*cos(Q1-Q2) + M*(L+LC)*cos(Q2);
D(3,3) = 2*M;

% D(1,1) = J + m*(L-LC2)^2;
% D(1,2) = -J - (L-LC2)^2*m + L*(L-LC2)*m*cos(Q1);
% D(2,1) = -J - (L-LC2)^2*m + L*(L-LC2)*m*cos(Q1);
% D(2,2) = 2*J + 2*L^2*m - 2*L*LC2*m + LC2^2*m + LC^2*M - 2*L*(L-LC2)*m*cos(Q1);

G = zeros(3,1);
G(1) = GRAV*M*(L-LC)*sin(Q1-Q2);
G(2) = -GRAV*M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
G(3) = 0;

% G(1) = GRAV*m*(L-LC2)*sin(Q1-Q2-gam);
% G(2) = -GRAV*(L-LC2)*m*sin(Q1-Q2-gam) - GRAV*(L*m+LC*M)*sin(Q2+gam);%-GRAV*M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));

C = zeros(3);

C(1,2) = -L*M*(L-LC)*sin(Q1)*U2;
C(2,1) = -L*M*(L-LC)*sin(Q1)*(U1-U2);
C(2,2) = L*M*(L-LC)*sin(Q1)*U1;
C(3,1) = -M*(L-LC)*sin(Q1-Q2)*U1 + 2*M*(L-LC)*sin(Q1-Q2)*U2;
C(3,2) = -M*(L-LC)*sin(Q1-Q2)*U2 - M*(L+LC)*sin(Q2)*U2;

% C(1,2) = L*(-L+LC2)*m*U2*sin(Q1); % -L*M*(L-LC)*sin(Q1)*U2;
% C(2,1) = -L*(L-LC2)*m*(U1 - U2)*sin(Q1) ;%-L*M*(L-LC)*sin(Q1)*(U1-U2);
% C(2,2) = L*(L-LC2)*m*U1*sin(Q1) ;%L*M*(L-LC)*sin(Q1)*U1;

B1 = [1; 0; 0];
B2 = [0; 0; 1];

% This might be for adding disturbance (TChen)
% J_com=zeros(2,2);
% J_com(1,1) = .5*(L-LC)*cos(Q1-Q2);
% J_com(1,2) = .5*(L+LC)*cos(Q2) - .5*(L-LC)*cos(Q1-Q2);
% J_com(2,1) = .5*(L-LC)*sin(Q1-Q2);
% J_com(2,2) = -.5*(L+LC)*sin(Q2) - .5*(L-LC)*sin(Q1-Q2);
% % J_com(2,1) = -.5*(L-LC)*sin(Q1-Q2);
% % J_com(2,2) = -.5*(L+LC)*sin(Q2) + .5*(L-LC)*sin(Q1-Q2);
%
% J_hip=zeros(2,2);
% J_hip(1,2) = L*cos(Q2);
% J_hip(2,2) = -L*sin(Q2);

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

% look at H at contact
dH = zeros(1,6);
dH(1,1) = 1;
dH(1,2) = -(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3);

% LfH = dH(1, 1:3) * qds;
% Below is equivalent, but save time
LfH = x(4)-(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3)*x(5);

dLfH = zeros(1,6);
dLfH(1,2) =  -(2*b2+6*b3*x(2)+12*b4*x(2)^2)*x(5);
dLfH(1,4) = 1;
dLfH(1,5) =  -(b1+2*b2*x(2)+3*b3*x(2)^2+4*b4*x(2)^3);
