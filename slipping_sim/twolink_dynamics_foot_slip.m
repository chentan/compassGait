function dx = twolink_dynamics_foot_slip(t,x,params)
% this is an extended state x = [q1 q2 xs q1d q2d xsd]
% dx = f(x) + g(x)*u + gx_fric*u_fric

%% Assume that the original HZD for no-slip situation is adopted
% ns: no slip
x_ns = [x(1:2); x(4:5)];
params.a;
[D,C,G,B,H,~,LfH,dLfH] = twolink_EOM(x_ns,params);

Fx_ns = D\(-C*x_ns(3:4)-G); 

Gx_ns = D\B; 

% Calculate Lie Derivatives needed for control
LgLfH_ns = dLfH*[zeros(2,1);Gx_ns];
Lf2H_ns = dLfH*[x_ns(3:4);Fx_ns];

% Calculate control input
Kd = 25;
Kp = 25; 
e = .2;
v = pd_control(Kp,Kd,e,H,LfH);

% Bernstein-Bhat controller (uses feedback linearization)
% Copied / Modded from Eric's 3-link
u = LgLfH_ns \ (v-Lf2H_ns);


%% Interaction with slippery ground
% dx = f(x) + g(x)*u + gx_fric*u_fric
coeff_of_fric = params.fric_coeff;

[D,C,G,B1,B2,~,~,~,~] = twolink_EOM_foot_slip(x,params);

Fx = D\(-C*x(4:6)-G); 

Gx = D\B1; 
Gx_fric = D\B2; % due to the friction

% Compute the Friction Force
% FN = Const+U1p_coef*u1p+U2p_coef*u2p;
[Const,U1p_coef,U2p_coef] = stance_force_two_link_foot_slip(x);
% Assume always in contact between the foot and the ground and FN >= 0
% Need to reconsider about this part
% u_fric = -(coeff_of_fric)*FN 

if x(6) == 0
    x_check = [x(1:2)',x(4:5)'];
    dx_check = twolink_dynamics(0,x_check',params);
    [FT_check,~] = stance_force_two_link(x_check,dx_check);
    fric_sign = sign(FT_check);
else
    fric_sign = -sign(x(6));
end

% u_fric = fric_sign * coeff_of_fric * (Const + U1p_coeff * U1p + U2p_coeff * U2p);
% u_fric = fric_sign * coeff_of_fric * Const + 
% [fric_sign * coeff_of_fric* U1p_coeff, fric_sign * coeff_of_fric* U2p_coeff, 0] * 
% [U1p; U2p; ddxs];

B = Fx + Gx * u + Gx_fric * (fric_sign * coeff_of_fric * Const);
A = eye(3) - Gx_fric * [fric_sign * coeff_of_fric* U1p_coef, ...
    fric_sign * coeff_of_fric* U2p_coef, 0];
dx(1:3) = x(4:6);
dx(4:6) = linsolve(A,B);

%% Reorient
dx=dx';