function [H,dH,LfH,dLfH] = mod_H_fun(x,params);
% Load params
a=params.a;
[L,LC,M,J,GRAV] = parameters();

%% do the conversion from s to theta space
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


% b0=.829025;
% b1=1.8554;
% b2=-20.8159;
% b3=2.87177;
% b4=87.1477;
%% Compute angle addition
% JN added display commands to see what happens at disturbance -
% start_angle updates from x0(2) to x(2)
%disp(params.control.start_angle)
%disp(params.control.end_angle)
%disp(params.control.add_angle)

%params.control.start_angle = 0.2244 % JN added to test effect of constant start angle
%params.control.add_angle = 10*pi/180
m = -params.control.add_angle/(params.control.end_angle - params.control.start_angle);
x(2); % JN added for debugging
y = m.*(x(2) - params.control.start_angle) + params.control.add_angle;


%% add it to the H, dH etc, remain unchanged...

%  H = [x(1) + (- y - b0-b1*x(2)-b2*x(2)^2-b3*x(2)^3-b4*x(2)^4);];
H = [x(1) + (- b0-b1*x(2)-b2*x(2)^2-b3*x(2)^3-b4*x(2)^4);];

 
 dH=zeros(1,4);
 dH(1,2) = - (b1 + 2*b2*x(2) + 3*b3*x(2)^2 + 4*b4*x(2)^3);
 dH(1,1) = 1;

 LfH = [x(3) - (b1 + 2*b2*x(2) + 3*b3*x(2)^2 + 4*b4*x(2)^3)*x(4)];
 
 dLfH = zeros(1,4);
 dLfH(1,2) =  -(2*b2 + 6*b3*x(2) + 12*b4*x(2)^2)*x(4);
 dLfH(1,4) =  -(b1 + 2*b2*x(2) + 3*b3*x(2)^2 + 4*b4*x(2)^3);
 dLfH(1,3) = 1; 