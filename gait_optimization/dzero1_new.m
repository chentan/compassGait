function [xi10,xi20,dzero,zeta_star_2,maxVzero] = dzero1_new(params)
a=params.a;

a0 = a(1);
a1 = a(2);
a2 = a(3);
a3 = a(4);
a4 = a(5);

[L,LC,M,J,GRAV] = parameters();

% From symmetry and etc.
THp=-a4/2;
THm=-THp;

H0=[1 0];
c = [0 1];
H=[H0;c];

q_plus = inv(H)*[a0;THp];
q_minus = inv(H)*[a4;THm];
% a0 = swing leg
% a4/2 = stance
 
% try xminus / plus from calc above, try velocities zero?
 xminus = [q_minus' 0 0]'; 
 xplus = [  q_plus' 0 0]'; % JN - x starts here (nope, in find_parameters). Add another term?
 
 % changed to params
[D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(xminus,params); % JN - x fed into EOM here, should add term to B?
    [x0,delta_dq,F2]=impact_2link(xminus);

Q1=xminus(1);
Q2=xminus(2);
gamma0m=[(-1).*J+(-1).*(L+(-1).*LC).*M.*(L+(-1).*LC+(-1).*L.*cos(Q1)),2.* ...
  J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+(-2).*L.*(L+(-1).*LC).*cos(Q1) ...
  )];
lambda_bar_dq = inv([dH(1:2);gamma0m])*[0;1];
Q1=x0(1); % xplus
Q2=x0(2); % xplus
gamma0p=[(-1).*J+(-1).*(L+(-1).*LC).*M.*(L+(-1).*LC+(-1).*L.*cos(Q1)),2.* ...
  J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+(-2).*L.*(L+(-1).*LC).*cos(Q1) ...
  )];
dzero=gamma0p*delta_dq*lambda_bar_dq;
d2zero=dzero^2;
% Checks with Eric's book.  Table 6.1, updated from errata.

% qp = -pi/14;
% qm = pi/14;

qp = THp;
qm = THm;
% Conversion to b0 + b1*x + b2*x^2 +b3*x^3 +b4*x^4
b = [(qm+(-1).*qp).^(-4).*(a0.*qm.^4+qp.*((-4).*a1.*qm.^3+qp.*(6.*a2.* ...
    qm.^2+(-4).*a3.*qm.*qp+a4.*qp.^2))),(-4).*(qm+(-1).*qp).^(-4).*( ...
    a0.*qm.^3+(-1).*a1.*qm.^2.*(qm+3.*qp)+qp.*(3.*a2.*qm.*(qm+qp)+qp.* ...
    (a4.*qp+(-1).*a3.*(3.*qm+qp)))),6.*(qm+(-1).*qp).^(-4).*(a0.* ...
    qm.^2+a2.*qm.^2+4.*a2.*qm.*qp+(-2).*a3.*qm.*qp+a2.*qp.^2+(-2).* ...
    a3.*qp.^2+a4.*qp.^2+(-2).*a1.*qm.*(qm+qp)),(-4).*(qm+(-1).*qp).^( ...
    -4).*(a0.*qm+3.*a2.*qm+(-1).*a3.*qm+3.*a2.*qp+(-3).*a3.*qp+a4.*qp+ ...
    (-1).*a1.*(3.*qm+qp)),(a0+(-4).*a1+6.*a2+(-4).*a3+a4).*(qm+(-1).* ...
    qp).^(-4)];
b0=b(1);
b1=b(2);
b2=b(3);
b3=b(4);
b4=b(5);

vzerofun =@(Q2)(-1).*GRAV.*M.*(2.*J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+(-2).*L.*( ...
  L+(-1).*LC).*cos(b0+Q2.*(b1+Q2.*(b2+Q2.*(b3+b4.*Q2)))))+(-1).*(( ...
  -1).*b1+(-2).*b2.*Q2+(-3).*b3.*Q2.^2+(-4).*b4.*Q2.^3).*((-1).*J+( ...
  L+(-1).*LC).*M.*((-1).*L+LC+L.*cos(b0+Q2.*(b1+Q2.*(b2+Q2.*(b3+b4.* ...
  Q2))))))).*((-1).*(L+LC).*sin(Q2)+((-1).*L+LC).*sin(b0+Q2.*((-1)+ ...
  b1+b2.*Q2+b3.*Q2.^2+b4.*Q2.^3)));

THp=xplus(2);
THm=xminus(2);

% Integral
VzeroTHm = -quad(vzerofun,THp,THm, 1.e-6);

TH=linspace(THp,THm);
% Changing below
% Vzero(1) = 0, since TH(1) = THp, otherwise, quad gives singularity error
Vzero(1) = 0;
for i=2:length(TH)
    Vzero(i) = -quad(vzerofun,THp,TH(i));
end
maxVzero = max(Vzero);
zeta_star_2 = - VzeroTHm / (1-dzero^2);
% check1 should be less than 0
check1 = d2zero/(1-d2zero)*VzeroTHm + maxVzero;


% Explicit Qd apriori?  p. 155
% angle is clockwise so + angular momentum left to right?
% ERROR was here, forgot the 2.  blerg
gamma_star = sqrt(2*zeta_star_2);
q_dot_star = lambda_bar_dq*gamma_star;

% Find ICs for xis
R = [[-1 0];[-1 1]];
xi10 = c*R*xminus(1:2);
% now find for xi2, see EBC2) p. 155
foo = [xminus(1:2);q_dot_star];
[foo_new,delta_dq,F2]=impact_2link(foo);
Q1=foo_new(1);
Q2=foo_new(2);
Q1d=foo_new(3);
Q2d=foo_new(4);

gamma = Q1d.*((-1).*J+(-1).*(L+(-1).*LC).*M.*(L+(-1).*LC+(-1).*L.*cos(Q1)) ...
  )+Q2d.*(2.*J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+(-2).*L.*(L+(-1).* ...
  LC).*cos(Q1)));

%Try gamma on Za...  Seems to have same result ...
% Of course they have the same result
gamma = (b1+2.*b2.*Q2+3.*b3.*Q2.^2+4.*b4.*Q2.^3).*Q2d.*((-1).*J+(-1).*(L+( ...
  -1).*LC).*M.*(L+(-1).*LC+(-1).*L.*cos(b0+b1.*Q2+b2.*Q2.^2+b3.* ...
  Q2.^3+b4.*Q2.^4)))+Q2d.*(2.*J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+( ...
  -2).*L.*(L+(-1).*LC).*cos(b0+b1.*Q2+b2.*Q2.^2+b3.*Q2.^3+b4.*Q2.^4) ...
  ));
% This doesn't seem to work
% Update:  10 Apr 2009
% It works when you use the right formula for gamma*... sigh
xi20 = gamma;