function optim_cot = fun(alpha)

params.vel.flag=0;
params.flag2=0;
params.flag3=0;

params.vel.dist = 0;
params.vel.cycle_percent = 0.16;

params.vel.dir = 2;  % What direction is disturbance aligned with (see v_disturbance)
params.flags.step_one = 0;    % Is it the first step?
params.flags.ca_control = 1;  % Is the ca controller on or off
params.flags.vel_dist = 0;    % Is there a velocity disturbance? 1 Yes 0 No
params.control.angle_threshold = pi/7; % JN I believe this is x(2) at touchdown;

% setup parameters
params.flags.step_one = 1 ;  %No need to run undisturbed step, JN changed 0 to 1
params.flags.use_control = 0; % Don't use any extra control


%% Optimization of gaits
stepnum = 1;
params.a_0 = [0 0 alpha 0.4488];
data = run_walker4(params.a_0,stepnum,params);

if data.ieout(end) == 2
%     disp('hip close to the ground');
    wspd = 0;
    cot = 1e10;
elseif data.ieout(end) == 3
%     disp('falling');
    wspd = 0;
    cot = 1e10;
elseif data.ieout(end) == 6
%     disp('falling');
    wspd = 0;
    cot = 1e10;
else
    [L,~,M,~,GRAV] = parameters();
    q1temp = data.xout(end,1);q2temp = data.xout(end,2);
    totlg = L*(sin(q2temp)+sin(q1temp-q2temp))+data.p_horz(end,1);
    wspd = totlg/data.tout(end);

    n = size(data.tout);
    W = 0;
    for i = 1:n-1
        W = W+abs(data.uout(i)*(data.xout(i+1,1)-data.xout(i,1)));
    end
    
    wkpermt = W/totlg;
    cot = wkpermt/M/GRAV;
end

fileID = fopen('exp.txt','a');
fprintf(fileID,'%12.5f %12.5f %12.5f %12.5f\n',alpha(1),alpha(2), ...
    wspd,cot);
fclose(fileID);

optim_cot = cot;