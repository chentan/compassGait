function [angM,dangM] = angMomentum(t,x,params)

% [L,LC1,LC2,M,m,J,GRAV,gam] = parameters();
% LC =LC1;
[L,LC,M,J,GRAV] = parameters();

Q1=x(1);
Q2=x(2);
Q1d = x(3);
Q2d = x(4);


[D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(x,params);
Fx = inv(D)*(-C*x(3:4)-G);
Gx = inv(D)*B;

LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];

 [v]=control_two_link(H,LfH);
 u = inv(LgLfH)*(v-Lf2H);

angM = (-1/2).*M.*(2.*L.*LC.*Q1d+(-1).*L.^2.*(Q1d+(-2).*Q2d)+(-1).* ...
  LC.^2.*(Q1d+(-2).*Q2d)+(L.^2+(-1).*LC.^2).*(Q1d+(-2).*Q2d).*cos( ...
  Q1));

dtemp = [(1/2).*(L.^2+(-1).*LC.^2).*M.*(Q1d+(-2).*Q2d).*sin(Q1),0,(-1/2).* ...
  M.*((-1).*L.^2+2.*L.*LC+(-1).*LC.^2+(L.^2+(-1).*LC.^2).*cos(Q1)),( ...
  -1/2).*M.*(2.*L.^2+2.*LC.^2+(-2).*(L.^2+(-1).*LC.^2).*cos(Q1))];

dangM = dtemp*[x(3:4);Fx]; + dtemp*[zeros(2,1);Gx]*u;