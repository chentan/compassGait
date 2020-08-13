function u = switching_control(x,params)

params.a;
[D,C,G,B,H,~,LfH,dLfH] = twolink_EOM(x,params);

Fx = D\(-C*x(3:4)-G);

Gx = D\B;

LgLfH = dLfH*[zeros(2,1);Gx];
Lf2H = dLfH*[x(3:4);Fx];

Kd = 25;
Kp = 25; 
e = .2;
v = pd_control(Kp,Kd,e,H,LfH);

u = inv(LgLfH)*(v-Lf2H);