function Q = zero2full(xi,a)
% Function takes the outs from zerodyn2
% and converts it to the full model...
% Q(:,2) = xi(:,1)
[L,LC,M,J,GRAV] = parameters();

Q2 = xi(:,1);

a0 = a(1);
a1 = a(2);
a2 = a(3);
a3 = a(4);
a4 = a(5);

% Changed to correspond to a, dcp 21 APR 2009
qp = -a4/2;
qm = -qp;
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


Q1 = -(- b0-b1.*Q2-b2.*Q2.^2-b3.*Q2.^3-b4.*Q2.^4);

% Qhyp = -1.122; % Tan added, or Qhyp = Q1+3*Q2

xi1 = xi(:,1);

K1 = (2.*J+LC.^2.*M+M.*(L.^2+(L+(-1).*LC).^2+(-2).*L.*(L+(-1).*LC).* ...
  cos(b0+xi1.*(b1+xi1.*(b2+xi1.*(b3+b4.*xi1)))))+(-1).*((-1).*b1+(-2).* ...
  b2.*xi1+(-3).*b3.*xi1.^2+(-4).*b4.*xi1.^3).*((-1).*J+(L+(-1).*LC).* ...
  M.*((-1).*L+LC+L.*cos(b0+xi1.*(b1+xi1.*(b2+xi1.*(b3+b4.*xi1))))))).^( ...
  -1);

xi2=xi(:,2);

xi1d = K1.*xi2;
Q2d = xi1d;

% Q1d is determined by Q2d from the Bezier Poly
Q1d = (b1 + 2.*b2*Q2 + 3.*b3.*Q2.^2 + 4.*b4.*Q2.^3).*Q2d;

Q=[Q1 Q2 Q1d Q2d];