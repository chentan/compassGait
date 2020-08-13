function [name] = run_walker4(a_0,steps,params)
% Runs the 2-link gait associated w/ polynomial
%
% Inputs
% a - vector of coefficents of desired tracking polynomial
% length(a) = 5
% note that a(1) and a(2) can be equal to zero because
% the function will assign the correct values based on 
% the remaining parameters
%
% Outputs
% name - a structure of the pertinent data
% name is formated to be the input argument to all the plot routines
%

a = find_parameters(a_0); % Input coefficients
params.a = a; % JN Save these so that we can revert back after modifying
params.asave = a;

% changed for params
[xi10,xi20,dzero,zeta_star_2,maxVzero] = dzero1_new(params);

% now convert to full model...
Q = zero2full([xi10,xi20],a); % a is the alpha in 6.7 to calculate b term
x0 = Q'; 

% Ang. Mom. Mod
[angM,dangM] = angMomentum(0,x0,params);
x0(5) = angM;

tstart=0;
tfinal=3*steps;  % estimate how long for # of steps

% Regular options
options = odeset('Events','twolink_events','Refine',4,'RelTol',10^-6,'AbsTol',10^-6);

%% Initialize output data matricies

params.control.start_angle = x0(2); % from line ~170
params.control.end_angle = params.a_0(5)/2;
params.control.add_angle = 10*pi/180;

forces = [];
tout=[0];
xout=[x0'];
teout=[];
xeout=[];
ieout=[];
errorout = []; % JN added
p_horz=[0];
p_horz_1=0;
uout = switching_control(x0,params);

% Set to record data on 1st step
params.flags.step_one = 0;

%% Code to run for input # of steps if vel disturbance is present
if params.flags.vel_dist == 1
    loops =steps+2;
else
    loops = steps;
end

%% Main Simulation Loop
for n=1:(loops)
    % n % uncomment to display the loop count. Note that this increments after a velocity disturbance
    if n >=3
        params.a = params.asave; % Revert to original coefficients if running modification code
    end
    % Run a step
    [t,x,te,xe,ie] = ode45(@twolink_dynamics,[tstart tfinal],x0,options,params);
    
    
    % Save the data
    nt = length(t);
    tout = [tout; t(2:nt)];
    xout = [xout;x(2:nt,:)];
    teout = [teout; te]; % Events at tstart are never reported.
    xeout = [xeout; xe];
    ieout = [ieout; ie];
    u=[];
    
    for ii = 1:nt
        u(ii,:) = switching_control(x(ii,:)',params);
    end
    uout = [uout; u(2:nt)];
    
    % Update the time history of stance foot position
    [pFoot1,pFoot2,pHip,pTorso] = limb_pos(x(end,:),0);
    p_horz=[p_horz; p_horz_1*ones(length(x)-1,1)];
    
    % Record CoM position on first step and store it
    for ts = 1:nt
    [~,~,~,~,pCOM] = limb_pos(xout(ts,:),p_horz(ts));
%     [pFoot1,pFoot2,pHip,pTorso,pCOM] = limb_pos(xout(ts,:),p_horz(ts));
    name.pCOM(ts) = pCOM(1);
    end
    
    % JN - not really sure what this next block does
    % Possibly use this for control
    if params.flags.step_one == 1; 
        for m =1:size(xout,1)
            [scoord,wcoord,swterm,sv,ds,Y1,Yp,vCOM]=control_authority(xout(m,:),params);
            params.vCOM(m,:) = vCOM;
            params.xout(m,:) = xout(m,:);
            params.scoord(m,:) = scoord;
        end
        params.flags.step_one = 0;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Hack to use CA control on next step
%           params.flags.ca_control = 1;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % Depending on the step end condition, proceed:
    % No end event, probably a timeout
    if isempty(ie)
        break;
    end
    
    if ie(end)==4; % Simulation stopped at disturbance configuration
        disp('disturbance'); %
        x0 = x(nt,:);
        %         [xd_new,vCOM_new,vCOM,params] = v_disturbance(x0,params); % Find
        %         post disturbance velocity
        % use ~ to ignore outputs in 2009b+
        
        % Find post disturbance velocity
        
        [xd_new,temp,temp,params] = v_disturbance(x0,params);
               
        % Set ICs to reflect this
        x0(3:4) = xd_new;
        
% JN - this next block is the Virtual Constraint Modification code. Uncommenting
% this will run an optimization for the new polynomial coefficients after a
% disturbance. To work, a target touchdown state must be specified within
% the VC_mod function. That function contains several targets that I used.

%         disp('Start optimization');
%         guess = [1,1,1,1,1];
% 
%         optionsopt = optimoptions(@fmincon,'Algorithm','sqp');
%         optionsopt.MaxFunEvals = 2000;
%         [c,dist] = fmincon(@(c)VC_mod(c,params,x0,tstart, tfinal),guess,[],[],[],[],[-10000;-10000;-10000;-10000;-10000],[],[],optionsopt);
%         
%         params.a = c
%         disp('Resolved');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Turn on the switching control if used
        if params.flags.use_control == 1;
            params.flags.ca_control = 1;  % Turn on CA control
        else
            params.flags.ca_control = 0;
        end
        % Various Recovery Methods to try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Switch to a new parameter set
        %         a_0 = [0 0 1.8700 1.2091 0.4488]; %168
        %         a_0 = [0 0 1.2775 0.6200 0.4488]; % 50
        %         a_0 = [0 0 1.4900 0.9031 0.4488] % 100
        %         a_0 = [0 0 1.7600    1.1352 0.4488]; % 150
        %         a = find_parameters(a_0);
        %         params.a = a;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate a modified trajectory for recovery
        params.control.start_angle = x0(2);
        params.control.end_angle = params.a_0(5)/2;
               params.control.add_angle = 10*pi/180; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Disturbance has happened, don't run it again
        params.flags.vel_dist =0;
        
        % These are additional cut-offs to the swing leg motion
    elseif ie(end)==5 || ie(end)==6
        params.flags.ca_control = 0;  % Turn off CA control
        params.flags.vel_dist = 0;  % only disturb once
        x0 = x(nt,:);
        %         [angM,dangM] = angMomentum(0,x0,params);
        % x0(5) = angM;
        n=n-1;  % don't count step
        % Regular end of step conditions
    else
        % Calculate new stance foot x-position
        p_horz_1=p_horz(end)+pFoot2(1);
        
        % Run the impact map
        [x0,delta_dq,F2]=impact_2link(x(nt,:)');
        % Display impact information
        %disp(['impact FT/FN = ', num2str(abs(F2(1))/F2(2)),' velocity = ' num2str(pFoot2(1)/(tout(end)-tstart))])
        %JN COMMENTED ^^^ out
        
        % Angular Momentum Mod
        [angM,dangM] = angMomentum(0,x0,params);
        x0(5) = angM;
        
        % JN - This block calculates the error (H) between the desired 
        % relationship between x(1) and x(2) and the actual relationship 
        [D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(xout(length(xout),:),params);
        % H
        xout(length(xout),:);
        if n==2
        name.xend = xout(length(xout),:);
        end
        
        %turn off extra control at end of step
        if params.flags.vel_dist ==0
            %           params.flags.use_control = 0;
            params.flags.ca_control = 0;
            %disp('CA OFF==================')
        end
        
    end
    
    
    
    tstart=t(nt);
    if tstart >= tfinal
        break
        % Falling ?Back?
    elseif ie(end)==2
        disp('falling!');
        break;
        % Falling ?Forward?
    elseif ie(end)==3
        disp('falling!');
        break;
    end
end

%% Calculate some output data

% Forces on the stance foot
forces = zeros(length(xout),2);
for j=1:length(xout)
    dxout = twolink_dynamics(0,xout(j,:)',params);
    [FT,FN] = stance_force_two_link(xout(j,:),dxout);
forces(j,:) = [FT FN];
end

% Control Authority quantities
[scoord,wcoord,swterm,sv,ds,Y1,Yp,vCOM]=control_authority(xout,params);

% Calculate errors in each step

for iii = 1:length(xout) % JN added
   [D,C,G,B,H,dH,LfH,dLfH] = twolink_EOM(xout(iii,:)',params);
   error(iii,1) = H; 
end


%% Put into output structure
    name.a=a;
    name.tout=tout;
    name.xout=xout;
    name.uout=uout;
    name.teout=teout;
    name.xeout=xeout;
    name.ieout=ieout;
    name.scoord=scoord;
    name.wcoord=wcoord;
    name.swterm=swterm;
    name.Y1=Y1;
    name.Yp=Yp;
    name.vCOM=vCOM;
    name.sv=sv;
    name.ds=ds;
    name.forces=forces;
    name.p_horz=p_horz;
    name.params = params;
    name.error = error; 
    
return;