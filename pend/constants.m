% smapling period
%global T_samp 

jitter = 1e-8;

% PENDULUM
% T_samp = 1/5;
% x1_min = -2;
% x1_max = 2;
% x2_min = -2;
% x2_max = 2;
% u_min = -4;
% u_max = 4;
% x0 = [pi/2; 0]; % NOT HERE!
% xs = [0; 0];
% us = 0;
% Q = 10*eye(nx);
% R = 0.5*eye(nu);
% P = 50*eye(nx);
% delta_bar = 0.1; 

% CSTR
x1_min = 1;
x1_max = 3;
x2_min = 0.5;
x2_max = 2;
u_min = 3;
u_max = 25; 

xs = [2.14; 1.09];
us = 14.19;
Q = diag([0.2 1]); 
R = 0.1;
P = [14.46 13.56; 13.56 62.22];
T_samp = 30;
delta_bar = 0.01; 

nx = 2;
nu = 1;

x_min = [x1_min; x2_min];
x_max = [x1_max; x2_max];

% kernel function

% some neat colors
AZZURRO = [0.8500 0.3250 0.0980];
RIPEORANGE = [0 0.4470 0.7410];