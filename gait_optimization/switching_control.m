function u = switching_control(x,params)
[D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(x,params);
% [D,C,G,B,H,~,LfH,dLfH] = twolink_EOM(x,params);
% Fx = inv(D)*(-C*x(3:4)-G);
Fx = D\(-C*x(3:4)-G); % faster than inv(D)*

% Gx = inv(D)*B;
Gx =     D\B; % faster than inv(D)*
    
LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];

flag = params.flags.ca_control;
%% Check Flags to Determine which controller to use
% if flag ==1 && angle <= angle_threshold
%     u = ca_control  % JN dig around for this, try to find where the value of ca_control is set
% elseif flag == 1 angle > angle_threshold
%     %turn off the CA controller
%     % somehow need to return this to parameters
%     % use HZD
%     [v]=control_two_link(H,LfH);
%     u = inv(LgLfH)*(v-Lf2H);
% else
%     [v]=control_two_link(H,LfH);
%     u = inv(LgLfH)*(v-Lf2H);

if flag ==1;
%     K1 = .05;
%     K2 = pi/2;
%     K3 = 5;
    
    % base off of uncontrolled velocity error
%     [scoord,wcoord,swterm,sv,ds,Y1,Yp,vCOM]=control_authority(x',params);
%     
%     vel_error = (scoord - interp1(params.xout(:,2),params.scoord,x(2)));
%     
% %     vel_error = (x(4) - interp1(params.xout(:,2),params.xout(:,4),x(2)));
%     if isnan(vel_error);  % happens if config that outside normal range
%         vel_error = .15; % guess something
%     end
% %      if  x(1) < K2 ; %abs(vel_error) > K1 &&
%         u=K3 * abs(vel_error);
%      else
%          u = 0.1;
%      end
         %         [v]=control_two_link(H,LfH);
%         u = inv(LgLfH)*(v-Lf2H);
%     end

% %%% S function approach
% [H,dH,LfH,dLfH,LgLfH,Lf2H] = CA_H_fun(x,params);
% Kd = 6;
% Kp = 8;
% e=1;
% v = pd_control(Kp,Kd,e,H,LfH);
% 
% u = inv(LgLfH)*(v-Lf2H);
% if u > 10
%     u = 10;
% elseif u < -10
%     u = -10;
% end
% if u == NaN
%     disp('woah!')
% end
% %%%%%%%%%%

%%%% Use Q2d-Q2d_desired
% [H,dH,LfH,LgH] = CA_H_Q2d(x,params);
% Kp = 2;
% v = -Kp*H;
% u = inv(LgH)*(v-LfH);


% %%% s function
% [H,dH,LfH,dLfH,LgLfH,Lf2H] = CA_H_fun(x,params);
% Kd = 2.5;
% Kp = 1.5;
% e=.5;
% v = pd_control(Kp,Kd,e,H,LfH);
% 
% u = inv(LgLfH)*(v-Lf2H);

% Use modified HZD H values
[H,dH,LfH,dLfH] = mod_H_fun(x,params);
Fx = D\(-C*x(3:4)-G); % faster than inv(D)*
Gx =     D\B; % faster than inv(D)*
LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];
    [v]=control_two_link(H,LfH);
    u = inv(LgLfH)*(v-Lf2H);        % figured out what in here changed when CA is on (and there is a disturbance)?
    % Before disturbance:
    % H = -0.0873
    % dH = [1.0000  -3.7096  0  0]
    % LfH = 1.1293e-07
    % dLfH = [0  10.5368  1.0000  -3.7096]
    % Fx = [-8.2442 ; 0.8847]
    % Gx = [23.8658  -0.3479]
    % LgLfH = 25.1563
    % Lf2H = -8.8640
    % v = 13.5963
    % u = 0.8928
    
    % After disturbance:
    % H = -0.1745
    % dH = [1.0000  -3.7096  0  0]
    % LfH = -18.0794
    % dLfH = [0  213.7870  1.0000  -3.7096]
    % Fx = [13.7085  -0.1191]
    % Gx = [23.8658  -0.3479]
    % LgLfH = 25.1563
    % Lf2H = 1.1101e+03
    % v = 340.8070
    % u = -30.5804
    
    % params.vel.dist = 4 % JN changed
    % params.vel.cycle_percent = 0.50
    % use_control & ca_control = 1
    
    % if params.vel.dist = 5, after the disturbance:
    % v = 413.2268
    % u = -51.3861
    
else
    [v]=control_two_link(H,LfH);
    u = inv(LgLfH)*(v-Lf2H);
    u = 1;
    v = 1;
end
