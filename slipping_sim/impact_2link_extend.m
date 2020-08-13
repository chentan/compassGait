function [x_new,delta_dq,F2,indNegFn,indNegVn]=impact_2link_extend(x,params)
coeff_of_fric = params.fric_coeff;
scale_coeff = params.scale_coeff;

%% Modified 
delta_dq = zeros(1,1); % just to make it complete

[L,LC,M,J,~] = parameters();

Q1 = x(1);
Q2 = x(2);

% Page 56 in Westvelt's book
%% Need to check De to see whether it is also suitable for slip case
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

%
E2(1,1) = L*cos(Q1-Q2);
E2(1,2) = L*(cos(Q2)-cos(Q1-Q2));
E2(1,3) = 1;
E2(1,4) = 0;
E2(2,1) = L*sin(Q1-Q2);
E2(2,2) = -L*(sin(Q2)+sin(Q1-Q2));
E2(2,3) = 0;
E2(2,4) = 1;

% derivative of the x position of touch-down foot end wrt the coordinates qe
der_x = E2(1, :);
% derivative of the y position of touch-down foot end wrt the coordinates qe
der_y = E2(2, :);

R = [[-1 0];[-1 1]];

%%
% Assume no slips
% eqn 3.20
left_mat = [De, -E2'; E2, zeros(2,2)];
qe_dot_bf = [x(4:6); 0]; % qe = [qj; pe]
right_vec = [De * qe_dot_bf; zeros(2,1)];

% touch-down foot end x-velocity before impact
vx_td_bf = dot(qe_dot_bf, der_x);
vy_td_bf = dot(qe_dot_bf, der_y);

temp = left_mat \ right_vec;
temp_ = zeros(1,6);

% F2 = [F_T; F_N]
F2 = temp(5:6);

if abs(F2(1)/F2(2)) <= (coeff_of_fric) * scale_coeff
    % no slips
%     disp('No slip at impact');
    params.flags.foot_slip = 0;
    if F2(2) < 0
        indNegFn = 1;
        indNegVn = [];
        x_new = [];
        disp('Foot taking off..');
        return
    else
        xd = temp(1:2);
        if temp(4) < 0
            indNegFn = 0;
            indNegVn = 1;
            x_new = [];
            disp('Warning: the stance leg interacts with the ground!')
            return
        end
    end
else
    % slips
    disp('Slips at impact')
    params.flags.foot_slip = 1;
    % recompute the impact map for slip case
    left_mat_1 = [De, -E2'; der_y, zeros(1,2); zeros(1,4), -1, (coeff_of_fric*scale_coeff)];
    temp_1 = left_mat_1 \ right_vec;
    sw_vel_x_1 = dot(temp_1(1:4), der_x); % touch-down foot x-vel after impact
%     if sign(sw_vel_x_1) == -sign(temp_1(5))
    if sign(vx_td_bf ) == -sign(temp_1(5))
        temp_ = temp_1;
        F2 = temp_1(5:6);
        if F2(2) < 0
            indNegFn = 1;
            indNegVn = [];
            x_new = [];
            disp('Foot taking off..');
            return
        else
            xd = [temp_1(1:2); sw_vel_x_1];
            if temp_1(4) < 0
                indNegFn = 0;
                indNegVn = 1;
                x_new = [];
                disp('Warning: the stance leg interacts with the ground!')
                return
            end
        end
    else
        left_mat_2 = [De, -E2'; der_y, zeros(1,2); zeros(1,4), 1, (coeff_of_fric*scale_coeff)];
        temp_2 = left_mat_2 \ right_vec;
        sw_vel_x_2 = dot(temp_2(1:4), der_x);
%         if sign(sw_vel_x_2) == -sign(temp_2(5))
        if sign(vx_td_bf ) == -sign(temp_2(5))
%             disp('Good, friction and velocity directions checked to be opposite.');
            temp_ = temp_2;
            F2 = temp_2(5:6);
            if F2(2) < 0
                indNegFn = 1;
                indNegVn = [];
                x_new = [];
                disp('Foot taking off..');
                return
            else
                xd = [temp_2(1:2); sw_vel_x_2];
                if temp_2(4) < 0
                    indNegFn = 0;
                    indNegVn = 1;
                    x_new = [];
                    disp('Attention: the stance leg interacts with the ground!')
                    return
                end
            end          
        else
            disp('No solution!') % this should never happen
            return
        end
    end
end

% swap the coordinates
if params.flags.foot_slip == 0
    x_new = zeros(4,1);
    x_new(1:2) = R * x(1:2);
    x_new(3:4) = R * xd;
else
    x_new = zeros(6,1);
    x_new(1:2) = R * x(1:2);
    x_new(4:5) = R * xd(1:2);
    x_new(6) = xd(3);
end
indNegFn = 0;
indNegVn = 0;

% Check the conservation of momentum using the CoM
syms q1 q2 xst yst 
xCOM = xst + (LC*sin(q2)+L*sin(q2)+(L-LC)*sin(q1-q2))/2;
yCOM = yst + (LC*cos(q2)+L*cos(q2)-(L-LC)*cos(q1-q2))/2;
xcom_der_symbol = [diff(xCOM, q1), diff(xCOM, q2), diff(xCOM, xst), diff(xCOM, yst)];
ycom_der_symbol = [diff(yCOM, q1), diff(yCOM, q2), diff(yCOM, xst), diff(yCOM, yst)];

xcom_der = eval(subs(xcom_der_symbol, [q1, q2, xst, yst], [Q1, Q2, 0 ,0]));
ycom_der = eval(subs(ycom_der_symbol, [q1, q2, xst, yst], [Q1, Q2, 0 ,0]));

vxcom_bf = dot(xcom_der, qe_dot_bf);
vycom_bf = dot(ycom_der, qe_dot_bf);

if params.flags.foot_slip == 0
    qe_dot_af = temp(1:4);
    vxcom_af = dot(xcom_der, qe_dot_af);
    vycom_af = dot(ycom_der, qe_dot_af);
%     disp(['x - impulse: ', num2str(temp(5)), ' = momentum difference ', num2str(2*M*(vxcom_af - vxcom_bf))]);
%     disp(['y - impulse: ', num2str(temp(6)), ' = momentum difference ', num2str(2*M*(vycom_af - vycom_bf))]);   
else
    qe_dot_af = temp_(1:4);
    vxcom_af = dot(xcom_der, qe_dot_af);
    vycom_af = dot(ycom_der, qe_dot_af);
%     disp(['x - impulse: ', num2str(temp_(5)), ' = momentum difference ', num2str(2*M*(vxcom_af - vxcom_bf))]);
%     disp(['y - impulse: ', num2str(temp_(6)), ' = momentum difference ', num2str(2*M*(vycom_af - vycom_bf))]); 
end

%% Original

% [L,LC,M,J,~] = parameters();
% 
% Q1 = x(1);
% Q2 = x(2);
% vvector = x(4:5); % velocity just before imparct
% x_new = zeros(4,1);
% 
% % Page 56 in Grizzle's book
% De(1,1) = J + M*(L-LC)^2;
% De(1,2) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
% De(1,3) = M*(L-LC)*cos(Q1-Q2);
% De(1,4) = M*(L-LC)*sin(Q1-Q2);
% De(2,1) = -J - M*(L-LC)*(L-LC-L*cos(Q1));
% De(2,2) = 2*J + M*LC^2 + M*(L^2+(L-LC)^2-2*L*(L-LC)*cos(Q1));
% De(2,3) = M*(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2));
% De(2,4) = -M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
% De(3,1) = M*(L-LC)*cos(Q1-Q2);
% De(3,2) = M*(L*cos(Q2)+LC*cos(Q2)-(L-LC)*cos(Q1-Q2));
% De(3,3) = 2*M;
% De(3,4) = 0;
% De(4,1) = M*(L-LC)*sin(Q1-Q2);
% De(4,2) = -M*(L*sin(Q2)+LC*sin(Q2)+(L-LC)*sin(Q1-Q2));
% De(4,3) = 0;
% De(4,4) = 2*M;
% E2(1,1) = L*cos(Q1-Q2);
% E2(1,2) = L*(cos(Q2)-cos(Q1-Q2));
% E2(1,3) = 1;
% E2(1,4) = 0;
% E2(2,1) = L*sin(Q1-Q2);
% E2(2,2) = -L*(sin(Q2)+sin(Q1-Q2));
% E2(2,3) = 0;
% E2(2,4) = 1;
% 
% R = [[-1 0];[-1 1]];
% 
% % Velocities pre-swap
% temp_vector = inv([De, -E2';E2 zeros(2)])*[De;zeros(2,4)];
% temp2 = inv([De, -E2';E2 zeros(2)])*[De*[vvector;0;0];0;0]; % Eqn 3.20 in book
% d_bar_dq = temp_vector(1:4,1:2);
% % check
% 
% F2 = temp2(5:6,:);
% 
% delta_F2 = -inv(E2*inv(De)*E2')*E2*[eye(2);zeros(2)];
% 
% % d_bar_dq_2 = inv(De)*E2'*delta_F2+[eye(2);zeros(2)];
% % d_bar_dq_2 - d_bar_dq
% impulse_F = E2'*F2;
% 
% delta_dq = [R zeros(2)]*d_bar_dq;
% % Swap the legs
% x_new(1:2) = R* x(1:2);
% % x_new2(3:4) = R* temp_vector(1:2,1:2)*vvector;
% x_test = 0;
% x_new(3:4) = R* temp2(1:2);
% x_new;
% 
% %% Check whether slip happens at impact
% % COM position for computing COM velocity 
% % xCOM = (Lc*sin(q2)+L*sin(q2)+(L-Lc)*sin(q1-q2))/2;
% % yCOM = (Lc*cos(q2)+L*cos(q2)-(L-Lc)*cos(q1-q2))/2;
% 
% % COM velocity
% Q1d = x(4); Q2d = x(5); xs = x(6);
% vCOMx = (LC*cos(Q2)*Q2d+L*cos(Q2)*Q2d+(L-LC)*cos(Q1-Q2)*(Q1d-Q2d))/2+xs;
% vCOMy = (-LC*sin(Q2)*Q2d-L*sin(Q2)*Q2d+(L-LC)*sin(Q1-Q2)*(Q1d-Q2d))/2;
% 
% Q1_new = x_new(1); Q2_new = x_new(2); Q1d_new = x_new(3); Q2d_new = x_new(4);
% vCOMx_new = (LC*cos(Q2_new)*Q2d_new+L*cos(Q2_new)*Q2d_new+ ...
%     (L-LC)*cos(Q1_new-Q2_new)*(Q1d_new-Q2d_new))/2;
% vCOMy_new = (-LC*sin(Q2_new)*Q2d_new-L*sin(Q2_new)*Q2d_new+ ...
%     (L-LC)*sin(Q1_new-Q2_new)*(Q1d_new-Q2d_new))/2;
% 
% fric_coeff = (coeff_of_fric);
% if fric_coeff < abs((vCOMx_new-vCOMx)/(vCOMy_new-vCOMy))
%     disp(['Slip at impact: ',num2str(abs((vCOMx_new-vCOMx)/(vCOMy_new-vCOMy)))]);
% else
%     disp(['No slip at impact: ', num2str(abs((vCOMx_new-vCOMx)/(vCOMy_new-vCOMy)))]);
% end