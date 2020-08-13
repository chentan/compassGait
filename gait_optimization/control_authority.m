function [scoord,wcoord,swterm,sv,ds,Y1,Yp,vCOM]=control_authority(x,params)

Q1= x(:,1);
Q2= x(:,2);
Q1d= x(:,3);
Q2d= x(:,4);

[L,LC,M,J,GRAV] = parameters();

% Convert relative angle to absolute angle, and switch angle symbols
theta1 = pi/2-Q2;
theta2 = 3*pi/2+Q1-Q2;
theta1d = -Q2d;
theta2d = Q1d-Q2d;

% scoord = (-1).*2.^(-1/2).*(J+(L.^2+(-1).*L.*LC+LC.^2).*M+(-1).*L.*(L+(-1).* ...
%   LC).*M.*cos(Q1)).^(-1/2).*(J.*(Q1d+(-2).*Q2d)+M.*(L.^2.*(Q1d+(-2) ...
%   .*Q2d)+LC.^2.*(Q1d+(-2).*Q2d)+2.*L.*LC.*((-1).*Q1d+Q2d))+(-1).*L.* ...
%   (L+(-1).*LC).*M.*(Q1d+(-2).*Q2d).*cos(Q1));

scoord = ((J+(L.^2+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2)).*theta1d+...
    (J+(L-LC).^2.*M+L.*(L-LC).*M.*cos(theta1-theta2)).*theta2d)./(sqrt(2).*sqrt(...
    J+(L.^2-L.*LC+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2)));

% wcoord = (1/2).*Q1d.*(J+(L.^2+(-1).*L.*LC+LC.^2).*M+(-1).*L.*(L+(-1).*LC).* ...
%   M.*cos(Q1)).^(-1).*(J.^2+2.*J.*(L.^2+(-1).*L.*LC+LC.^2).*M+(L+(-1) ...
%   .*LC).^2.*(L.^2+LC.^2).*M.^2+(-1).*L.^2.*(L+(-1).*LC).^2.*M.^2.* ...
%   cos(Q1).^2);

wcoord = (sqrt((J+(L.^2-L.*LC+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2))./(...
    J.^2+2.*J.*(L.^2-L.*LC+LC.^2).*M+(L-LC).^2.*(L.^2+LC.^2).*M.^2-L.^2.*(L-LC).^2.*M.^2.*(cos(theta1-theta2)).^2))...
    .*(2.*J.^2+4.*J.*(L.^2-L.*LC+LC.^2).*M+(L-LC).^2.*(L.^2+2.*LC.^2).*M.^2-L.^2.*(L-LC).^2.*M.^2.*cos(2.*(theta1-theta2)))...
    .*(theta1d-theta2d))./(2.*sqrt(2).*(J+(L.^2-L.*LC+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2)));

% swterm = (-2).*L.*(L+(-1).*LC).*M.*((-2).*J.^2+(-4).*J.*(L.^2+(-1).*L.*LC+ ...
%   LC.^2).*M+(-1).*(L+(-1).*LC).^2.*(L.^2+2.*LC.^2).*M.^2+L.^2.*(L+( ...
%   -1).*LC).^2.*M.^2.*cos(2.*Q1)).^(-1).*sin(Q1);

swterm = -(L.*(L-LC).*M.*sqrt((J+(L.^2-L.*LC+LC.^2).*M+L.*(L-LC).*M.*...
    cos(theta1-theta2))./(J.^2+2.*J.*(L.^2-L.*LC+LC.^2).*M+(L-LC).^2.*(L.^2+LC.^2).*M.^2-...
    L.^2.*(L-LC).^2.*M.^2.*(cos(theta1-theta2)).^2)).*sin(theta1-theta2))./(sqrt(2).*(J+(L.^2-...
    L.*LC+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2)));

% sv = (-1).*2.^(-1/2).*GRAV.*M.*(J+(L.^2+(-1).*L.*LC+LC.^2).*M+(-1).*L.* ...
%     (L+(-1).*LC).*M.*cos(Q1)).^(-1/2).*((L+(-1).*LC).*sin(Q1+(-1).*Q2) ...
%     +(L+LC).*sin(Q2));

sv = (GRAV.*M.*((L+LC).*cos(theta1)+(L-LC).*cos(theta2)))./(sqrt(2).*sqrt(...
        J+(L.^2-L.*LC+LC.^2).*M+L.*(L-LC).*M.*cos(theta1-theta2)));

for i=1:size(x,1)
    [D,C,G,B,H,dH,LfH,dLfH,J_com,J_hip] = twolink_EOM(x(i,:),params);
    J_hip=J_com;
    Y1_temp = (inv(D)*[1;0])';
    Y1(i,:) = (J_hip*Y1_temp')'/norm(J_hip*Y1_temp');
%     Y1(i,:) = Y1_temp/norm(Y1_temp);
    Yp_temp = ([0;1]/sqrt([0,1]*D*[0;1]))';
    Yp(i,:) = (J_hip*Yp_temp')'/norm(J_hip*Yp_temp');
    %
    % find COM Velocity
    %
    vCOM_temp = J_com*[Q1d(i);Q2d(i)];
    vCOM(i,:) = [vCOM_temp',norm(vCOM_temp)];
%     Yp(i,:) = Yp_temp/norm(Yp_temp);
end

ds= -scoord.*wcoord.*swterm - sv;