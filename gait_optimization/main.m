% AN OPTIMIZATION PROGRAM TO FIND THE OPTIMAL GAIT CORRESPONDING TO A
% SPECIFIED SPEED

clear all; clc

% a file to record the optimization history
fileID = fopen('exp.txt','w');
fprintf(fileID,'%12s %12s %12s %12s\n','a2','a3','speed','cot');
fclose(fileID);

% gait parameters a2 and a3. Give it an initial guess
alpha_0 = [1,0.38];

ub = [2,2]; % upper bound
lb = [0,0]; % lower bound

options = optimset('Display','iter','Algorithm','active-set', 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
alpha = fmincon(@(alpha)fun(alpha),alpha_0,[],[],[],[],lb,ub,...
    @(alpha)nonlcon(alpha),options);