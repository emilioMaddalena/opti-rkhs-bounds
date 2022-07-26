% smapling period
%global T_samp 
T_samp = 1/5;

jitter = 1e-8;

% system constraints
x1_min = -3;
x1_max = 3;
x2_min = -2;
x2_max = 2;
u_min = -1;
u_max = 1;

nx = 2;
nu = 1;

x_min = [x1_min; x2_min];
x_max = [x1_max; x2_max];

% kernel function



