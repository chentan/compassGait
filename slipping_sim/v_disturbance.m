function [xd_new,vCOM_new,vCOM,params] = v_disturbance(x0,params)
dir = params.vel.dir;
dist = params.vel.dist;

[D,C,G,B,H,dH,LfH,dLfH,J_com,J_hip] = twolink_EOM(x0,params);
vCOM = J_com*x0(3:4)';

% the 2-norm of Y1 and Yp are 1
if dir == 1
    vCOM_new = vCOM+[dist;0];
elseif dir == 2
    vCOM_new = vCOM+[0;dist];
elseif dir == 3
    vCOM_new = vCOM+dist*vCOM/norm(vCOM);
end

xd_new = inv(J_com)*vCOM_new;
xd_new = xd_new';
