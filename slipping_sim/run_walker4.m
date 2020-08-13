function [name] = run_walker4(a_0, steps, params)

coeff_of_fric = params.fric_coeff;
scale_coeff = params.scale_coeff;

a = find_parameters(a_0);
params.a = a;
params.asave = a;

% changed for params
[xi10,xi20,dzero,zeta_star_2,maxVzero] = dzero1_new(params);

% initialization
name.comment = [];

% now convert to full model
Q = zero2full([xi10,xi20],a); % a is the alpha in 6.7 to calculate b term
x0 = Q';

tstart = 0;
tfinal = 5*steps;  % estimate how long for # of steps

% ode options
options = odeset('Events','twolink_events','Refine',4,'RelTol',10^-12,'AbsTol',10^-12);
options_foot_slip = odeset('Events','twolink_events_foot_slip','Refine',4,'RelTol',10^-12,'AbsTol',10^-12);

%% Initialize output data matricies
F2_out = []; % impact forces
forces_out = [0,0]; % forces during swing phase
tout = [0];
x0_trans = x0';
xout = [x0_trans(1,1:2),0,x0_trans(1,3:4),0]; % all saved xout have 6 dimensions
teout = [];
xeout = [];
ieout = [];
p_horz = [0];
p_horz_1 = 0; % initial horizontal foot position is 0
pCOM_out = [];
uout = [0];

%% Code to change input # of steps if vel disturbance is present
if params.flags.vel_dist == 1
    loops = steps + 2;
else
    loops = steps;
end

%% Main Simulation Loop
n = 1; % initialize the first step
while n <= loops
    
    % at each beginning, check whether to slip or not
    if params.flags.foot_slip == 0
        dx_init = twolink_dynamics(0,x0,params);
        [FT_init,FN_init] = stance_force_two_link(x0,dx_init);
        if abs(FT_init/FN_init) > coeff_of_fric * scale_coeff
            params.flags.foot_slip = 1;
            x0_temp = x0;
            % extend a 4-dimensional vector to a 6-dimensional vector
            x0 = [x0_temp(1:2,1)',p_horz_1,x0_temp(3:4,1)',0]';
%             disp(['Slip at ',num2str(tout(end)), 's'])
        end
    end
    
    % continue
    if params.flags.foot_slip == 0
        % No slip
        [t,x,te,xe,ie] = ode15s(@twolink_dynamics,[tstart tfinal],x0,options,params);
        
        % Save the data
        nt = length(t);
        tout = [tout;t(2:nt)];
        x_ext = [x(:,1:2),p_horz_1*ones(nt, 1),x(:,3:4),zeros(nt, 1)];
        xout = [xout; x_ext(2:nt,:)];
        teout = [teout; te];
        % xe could be empty when there is a timeout
        if isempty(xe)
            xe_ext = [];
        else
            [row_num_xe, ~] = size(xe);
            xe_ext = [xe(:,1:2), p_horz_1*ones(row_num_xe, 1), xe(:,3:4), ...
            zeros(row_num_xe, 1)];
        end
        xeout = [xeout; xe_ext];
        ieout = [ieout; ie];
        
        % Calculate control inputs
        u = zeros(nt,1);
        for ii = 1:nt
            u(ii,:) = switching_control(x(ii,:)',params);
        end
        
    else
        % Slip: initial condition always contains 6 states
        [t,x,te,xe,ie] = ode15s(@twolink_dynamics_foot_slip,[tstart tfinal],x0,options_foot_slip,params);
        
        nt = length(t);
        tout = [tout; t(2:nt)];
        xout = [xout; x(2:nt,:)];
        teout = [teout; te];
        xeout = [xeout; xe];
        ieout = [ieout; ie];
        
        % Calculate control inputs
        u = zeros(nt,1);
        for ii = 1:nt
            u(ii,:) = switching_control_foot_slip(x(ii,:)',params);
        end
    end
    
    % Save the control
    uout = [uout;u(1:nt-1)];
    
    % Horizontal position of the stance foot
    % can be affected by 1) slip 2) impact and switch stance foot
    p_horz_1 = xout(end,3);
    
    % Calculate stance foot position and save the history in an array
    [~,pFoot2,~,~,~] = limb_pos(x(end,:),0);
    if params.flags.foot_slip == 0
        % No slip
        p_horz = [p_horz;p_horz_1*ones(length(x)-1,1)];
    else
        % Slip
        p_horz = [p_horz;x(2:end,3)];
    end
    
    % Calculate forces and save it
    forces = zeros(nt,2);
    if params.flags.foot_slip == 0
        for ii = 1:nt
            dxout = twolink_dynamics(0,x(ii,:)',params);
            [forces(ii,1),forces(ii,2)] = stance_force_two_link(x(ii,:),dxout);
        end
        forces_out = [forces_out;forces(1:end-1,:)];
        if ismember(1,forces(:,2)<0)
            name.comment = 'NegativeFn'; % infeasible for negative support forces
            break
        end
    else
        for ii = 1:nt
            dxout = twolink_dynamics_foot_slip(0,x(ii,:)',params);
            FN = stance_force_two_link_foot_slip_after(x(ii,:),dxout);
            % Determine both the direction and magnitude of friction forces
            forces(ii,2) = FN;
            if x(ii,6) == 0
                x_check = [x(ii,1:2),x(ii,4:5)];
                dx_check = twolink_dynamics(0,x_check',params);
                [FT_check,~] = stance_force_two_link(x_check,dx_check);
                fric_sign = sign(FT_check);
                forces(ii,1) = fric_sign*coeff_of_fric*FN;
            else
                forces(ii,1) = -sign(x(ii,6))*coeff_of_fric*FN;
            end
        end
        forces_out = [forces_out;forces(1:end-1,:)];
        if ismember(1,forces(:,2)<0)
            name.comment = 'NegativeFn';
            break
        end
    end
       
    % Depending on the step end condition, proceed:
    % No end event, probably a timeout
    if isempty(ie)
        name.comment = 'TimeOut-NoEnd';
        break;
    end
    
    tstart=t(nt);
    if tstart >= tfinal % swing
        name.comment = 'TimeOut-Swing';
        break
    elseif ie(end) == 2
        name.comment = 'Failure-ClosetoGnd';
        %                 disp('Too close to the ground!');
        break;
    elseif ie(end) == 3
        name.comment = 'Failure-FallBack';
        %                 disp('falling backwards!');
        break;
    elseif ie(end) == 6
        name.comment = 'Failure-LargeVel';
        %                 disp('integration issue! too large joint velocity');
        break;
    end
    
    name.comment = 'Success';
    % Disturbance
    % Assume disturbance occurs when walking with no slip and occurs for once
    if ie(end) == 4 % Simulation stopped at disturbance configuration
        disp('disturbance'); %
        x0 = x(nt,:);
        % Find post disturbance velocity
        [xd_new,~,~,params] = v_disturbance(x0,params);
        % Set new ICs
        x0(3:4) = xd_new;
        x0 = x0';
        params.flags.vel_dist = 0;
        
    % Foot Slip
    elseif ie(end)== 5
        if params.flags.foot_slip == 0
            % Switch to foot slip model
%             disp(['Slip starts at ',num2str(t(end)),' s'])
            params.flags.foot_slip = 1;
            % Set ICs
            x0 = [x(nt,1:2),p_horz_1,x(nt,3:4),0]';
        else
            % Switch back to no-slip model
%             disp(['Slip ends at ',num2str(t(end)),' s'])
            params.flags.foot_slip = 0;
            % set ICs
            x0 = [x(nt,1:2),x(nt,4:5)]';
        end
        
        % Regular end of step conditions
    else
        % After impact, the second step simulates on ground with normal
        % coefficient of friction
        if params.vary_fric == 1
            coeff_of_fric = 1; % for the current function zone
            params.fric_coeff = 1;
            params.vary_fric = 0;
        end
        
        % Calculate new stance foot x-position
        p_horz_1 = p_horz_1+pFoot2(1); % switch stance foot
        
        % Run the impact map
        [x0,~,F2,indNegFn,indNegVn] = impact_2link_extend(xout(end,:)',params);
        n = n + 1;
        if indNegFn == 1
            name.comment = 'NegativeFn';
            break
        elseif indNegVn == 1
            name.comment = 'NegativeVn';
            break
        else
            F2_out = [F2_out; F2'];
            if size(x0,1) == 6
                % slip occurs at impact
                params.flags.foot_slip = 1;
                % x0 is 6 by 1 vector
                x0 = [x0(1:2);p_horz_1;x0(4:6)];
            else
                % no slip occurs at impact
                % x0 is  4 by 1 vector
                params.flags.foot_slip = 0;
            end
        end
    end
end

% Calculate CoM position and save it
nt = length(tout);
pCOM_out = zeros(nt,2);
for ts = 1:nt
    [~,~,~,~,temp] = limb_pos(xout(ts,:),p_horz(ts));
    pCOM_out(ts,:) = temp';
end

% Put into output structure
name.a=a;
name.tout=tout;
name.xout=xout;
name.uout=uout;
name.teout=teout;
name.xeout=xeout;
name.ieout=ieout;
name.forces_out=forces_out; % contact forces
name.p_horz=p_horz;
name.pCOM_out = pCOM_out;
name.params = params;
name.F2_out = F2_out; % implulse at impact

return;