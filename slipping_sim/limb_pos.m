function [pFoot1,pFoot2,pHip,pTorso,pCOM] = limb_pos(x,ph)

[L,LC,M,J,GRAV] = parameters();

Q1 = x(1);
Q2 = x(2);

pFoot1  = [ph;0];  %  Check westervelt's code about this
pFoot2  = pFoot1 + [L*(sin(Q2)+sin(Q1-Q2)); L*(cos(Q2)-cos(Q1-Q2))];
pHip    = pFoot1 + [L*sin(Q2); L*cos(Q2)];
pTorso = pHip;
pCM1 = pFoot1 + LC*[sin(Q2); cos(Q2)];
pCM2 = pHip + (L-LC)*[sin(Q1-Q2); -cos(Q1-Q2)];
pCOM = 1/(2*M)*(pCM1*M + pCM2*M);
