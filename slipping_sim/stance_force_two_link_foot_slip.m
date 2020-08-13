function [Const,U1p_coef,U2p_coef] = stance_force_two_link_foot_slip(x)
[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);
U1 = x(4);
U2 = x(5);

% FT = -M*(L*sin(Q2)*U2^2+LC*sin(Q2)*U2^2+(L-LC)*sin(Q1-Q2)*(U1-U2)^2-(L-LC)*cos(Q1-Q2)*U1p-(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2))*U2p);
% FN = M*(2*GRAV+(L-LC)*cos(Q1-Q2)*(U1-U2)^2+(L-LC)*sin(Q1-Q2)*U1p-L*cos(Q2)*U2^2-LC*cos(Q2)*U2^2-(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2))*U2p);

Const = M*2*GRAV+M*(L-LC)*cos(Q1-Q2)*(U1-U2)^2-M*L*cos(Q2)*U2^2-M*LC*cos(Q2)*U2^2;
U1p_coef = M*(L-LC)*sin(Q1-Q2);
U2p_coef = -M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
