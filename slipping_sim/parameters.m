function [L,LC,M,J,GRAV] = parameters(walker)
if nargin == 0
L       =     1.0;    %   leg length, m
LC      =     0.8;    %   Leg CoM Location, m
M       =     5;    %   Leg mass, kg
J       =     0.6;   %   Leg Interia about CoM, kg-m^2 JN: changed from 0.03   
GRAV    =     9.81;   %   accel due to grav, m/s^2

return;
end
switch lower(walker)
    case{'simple','wisse'}
L       =     1.0;    %   leg length, m
LC      =     1.0;    %   CoM Location, m
M       =     1.0;    %   Leg mass, kg
J       =     0.00;   %   Leg Interia about CoM, kg-m^2   
GRAV    =     1.00;   %   accel due to grav, m/s^2
end

