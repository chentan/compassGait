function [x_new,delta_dq,F2] = impact_2link(x,params)
% [L,LC1,LC2,M,m,J,GRAV,gam] = parameters();
% % Hack for uneven legs
% LC = LC1; 
[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);
vvector = x(3:4); % velocity just before imparct
x_new = zeros(4,1);

% Page 56 in Grizzle's book
% Grizzle they used the stance foot position as the pe
De(1,1) = J + M*(L-LC)^2;
De(1,2) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
De(1,3) = M*(L-LC)*cos(Q1-Q2);
De(1,4) = M*(L-LC)*sin(Q1-Q2);
De(2,1) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
De(2,2) = 2*J + M*LC^2 + M*(L^2+(L-LC)^2-2*L*(L-LC)*cos(Q1));
De(2,3) = M*(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2));
De(2,4) = -M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
De(3,1) = M*(L-LC)*cos(Q1-Q2);
De(3,2) = M*(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2));
De(3,3) = 2*M;
De(3,4) = 0;
De(4,1) = M*(L-LC)*sin(Q1-Q2);
De(4,2) = -M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
De(4,3) = 0;
De(4,4) = 2*M;
E2(1,1) = L*cos(Q1-Q2);
E2(1,2) = L*(cos(Q2)-cos(Q1-Q2));
E2(1,3) = 1;
E2(1,4) = 0;
E2(2,1) = L*sin(Q1-Q2);
E2(2,2) = -L*(sin(Q2)+sin(Q1-Q2));
E2(2,3) = 0;
E2(2,4) = 1;

R = [[-1 0];[-1 1]];

% Velocities pre-swap
temp_vector = inv([De, -E2';E2 zeros(2)])*[De;zeros(2,4)];
temp2 = inv([De, -E2';E2 zeros(2)])*[De*[vvector;0;0];0;0]; % Eqn 3.20 in book
d_bar_dq = temp_vector(1:4,1:2);
% check

F2 = temp2(5:6,:);

% Eqn 3.23
delta_F2 = -inv(E2*inv(De)*E2')*E2*[eye(2);zeros(2)];

% Eqn 3.24
d_bar_dq_2 = inv(De)*E2'*delta_F2+[eye(2);zeros(2)];

% d_bar_dq_2 - d_bar_dq
impulse_F = E2'*F2;

delta_dq = [R zeros(2)]*d_bar_dq;
% Swap the legs
x_new(1:2) = R*x(1:2);
% x_new2(3:4) = R* temp_vector(1:2,1:2)*vvector;
x_test = 0;
x_new(3:4) = R*temp2(1:2);
x_new;