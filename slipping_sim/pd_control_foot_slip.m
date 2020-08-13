function v = pd_control_foot_slip(Kp,Kd,e,H,LfH)
% Kp = propotional gain
% Kd = derivative gain
% e = tuning parameter
% v = -(1/e*Kd*LfH + 1/e^2*Kp*H);
v = -(1/e*Kd*LfH + 1/e^2*Kp*H);