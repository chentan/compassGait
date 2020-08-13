function [FT,FN] = stance_force_two_link(x,dx)
[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);
U1 = x(3);
U2 = x(4);
U1p = dx(3);
U2p = dx(4);

FT = -M*(L*sin(Q2)*U2^2+LC*sin(Q2)*U2^2+(L-LC)*sin(Q1-Q2)*(U1-U2)^2-(L-LC)*cos(Q1-Q2)*U1p-(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2))*U2p);
FN = M*(2*GRAV+(L-LC)*cos(Q1-Q2)*(U1-U2)^2+(L-LC)*sin(Q1-Q2)*U1p-L*cos(Q2)*U2^2-LC*cos(Q2)*U2^2-(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2))*U2p);