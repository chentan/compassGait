function u = switching_control_foot_slip(x,params)

coeff_of_fric = params.fric_coeff;

[D,C,G,B1,B2,H,dH,LfH,dLfH] = twolink_EOM_foot_slip(x,params);

Fx = D\(-C*x(4:6)-G); 

Gx = D\B1; 
Gx_fric = D\B2; % due to the friction

% Calculate Lie Derivatives needed for control
LgLfH = dLfH*[zeros(3,1);Gx];
LgLfH_fric = dLfH*[zeros(3,1);Gx_fric];

% Compute the Friction Force
[Const,U1p_coef,U2p_coef] = stance_force_two_link_foot_slip(x);

Lf2H = dLfH*[x(4:6);Fx];

%% Calculate control input
Kd = 25;
Kp = 25; 
e = .2;
v = pd_control_foot_slip(Kp,Kd,e,H,LfH);

if x(6) == 0
    x_check = [x(1:2)',x(4:5)'];
    dx_check = twolink_dynamics(0,x_check',params);
    [FT_check,~] = stance_force_two_link(x_check,dx_check);
    fric_sign = sign(FT_check);
else
    fric_sign = -sign(x(6));
end

%% Compute dx
dx(1:3) = x(4:6);

Const1 = Fx;
Const2 = Gx*(inv(LgLfH)*(v-Lf2H)-inv(LgLfH)*LgLfH_fric*(fric_sign*(coeff_of_fric)*Const));
Const3 = Gx_fric*(fric_sign*(coeff_of_fric)*Const);

Coef_B = Gx*[-inv(LgLfH)*LgLfH_fric*(fric_sign*(coeff_of_fric)*U1p_coef), ...
    -inv(LgLfH)*LgLfH_fric*(fric_sign*(coeff_of_fric)*U2p_coef),0];
Coef_C = Gx_fric*[fric_sign*(coeff_of_fric)*U1p_coef, ...
    fric_sign*(coeff_of_fric)*U2p_coef,0];

A = eye(3)-Coef_B-Coef_C;
B = Const1+Const2+Const3;
dx(4:6) = linsolve(A,B);

%% Calculate the control input
u1p = dx(4); u2p = dx(5);
FN = Const+U1p_coef*u1p+U2p_coef*u2p;
u_fric = fric_sign*(coeff_of_fric)*FN;
u = inv(LgLfH)*(v-Lf2H-LgLfH_fric*u_fric);