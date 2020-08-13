% test pick parameters
function a_new = find_parameters(a)
%
% Functions finds a0 and a1 from a3 a4
% Input 5 vector of a values, a0 a1 will be over-written
%
% Function from modifiying / clean up of test_pick_parameters
%
%David C Post
% 
%
a0 = a(1);
a1 = a(2);
a2 = a(3);
a3 = a(4);
a4 = a(5);

H0=[1 0];
c = [0 1];
H=[H0;c];
q2 = a4/2;
q1 = a4;

dq1 = 0;
dq2 = 0;
x=[q1,q2,dq1,dq2]';

% a4=q1;
th_minus=q2;

M=4;
% a3=0.8
% Try different position to check J
% a3=1.1;
% [a0,a1,th_plus]=corollary61(H,a4,th_minus,x)
% function [a0,a1,th_plus]=corollary61(H,a4,th_minus,x)
R = [[-1 0];[-1 1]];
temp=H*R*inv(H)*[a4;th_minus];
a0=temp(1);
th_plus=temp(2);
[x_new,delta_dq,F2]=impact_2link(x);
wa_minus = inv(H)*[M/(th_minus - th_plus)*(a4-a3);1];
dq2 = 1;

a1 = H0*delta_dq*wa_minus*(th_minus - th_plus)/M*inv(c*delta_dq*wa_minus)+ a0;

a_new=[a0 a1 a2 a3 a4];