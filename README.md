This repository contains two folders: slipping_sim and gait_optimization. All codes run in Matlab.

# Folder slipping_sim 
1. It provides the code to animate a feasible gait (obtained from no-slip rough surfaces) walking on slippery ground.

2. Run the main.m function to start the program. In the main.m function, the value of params.fric_coeff can be changed to vary the ground friction. 

3. Output: When the command window prints out "Success", it represents that the feasible gait is success (with or without slipping). When it prints out "Fall backward", it represents that the feasible gait fails on slippery ground by falling backward (violating foward motion constraint). When it prints out "Negative force required", it represents that the feasible gait fails by requiring negative contactforce (violating no-flight-phase constraint). An animation will start after printing out the message.

# Folder gait_optimization
1. It provides the code to obtain the most energetically-efficient gait with a speed specified.

2. Run the main.m function in the folder to start the program. 

3. Output: Matlab will display the solution process when it is searching for the optimal gait. A file exp.txt will be generated after running the program. This file record the history of the gait parameters, speed, and cost of transport (CoT) during the optimization process.
 
Author contact information:
Tan Chen tchen8@nd.edu
