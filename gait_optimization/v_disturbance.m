function [xd_new,vCOM_new,vCOM,params] = v_disturbance(x0,params)
dir = params.vel.dir;

[scoord,wcoord,swterm,sv,ds,Y1,Yp,vCOM]=control_authority(x0,params);
[D,C,G,B,H,dH,LfH,dLfH,J_com,J_hip] = twolink_EOM(x0,params);

% the 2-norm of Y1 and Yp are 1
if dir == 1
    dist_dir = Y1;
elseif dir==2
    dist_dir = Yp;
elseif dir==3;
    dist_dir = (Y1+Yp)/norm(Y1+Yp);
end

vCOM_temp = vCOM(1:2) + params.vel.dist*dist_dir;
vCOM_new = [vCOM_temp, norm(vCOM_temp)];
xd_new = inv(J_com)*vCOM_new(1:2)';
params.vel.B2=swterm;
params.vel.G=-sv;
