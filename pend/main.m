clear all
close all
clc

addpath('./utils')

rng(1)
global T_samp 

constants

kernel = @(x1,x2,lengthscale) exp(-dist(x1,x2').^2 / (2*lengthscale^2));
dynamics = @(t,x,u) cstr(t,x,u);

% PART1

% extracting a dataset from the system
N = 8;
D = 400*ones(1,N);
connected = false; 
visuals = false;
dataset = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, 0, 'grid');

% estimate the kernel lengthscales
lengthscale_range = 2:0.1:5;
folds = 5;
lambda = 0.0000001;
[lengthscales, losses] = cv(folds, kernel, lambda, dataset, lengthscale_range);

% estimate ground-truths RKHS norms
gammas = estimate_rkhs(dataset, @(x1,x2,lengthscale) kernel(x1,x2,lengthscale), lengthscales);

gammas
lengthscales
%%
% PART2
% Learning a 1-step ahead KRR model for the system

% extracting a dataset from the system
connected = false; 
visuals = false;    
dataset = collect_data(400, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'grid');

for nx = 1:nx
    
    X{nx} = dataset{nx,1}(:,1:end-1);
    y{nx} = dataset{nx,1}(:,end);

    jitter = 1e-8;
    n_data = size(X{nx},1);
    K{nx} = kernel(X{nx},X{nx},lengthscales(nx,1)) + jitter*eye(n_data);
    alpha{nx} = (K{nx} + n_data*lambda*eye(n_data))\y{nx}; 
    
    kernel_casadi = @(x1,x2,lengthscale) exp(-diag((x1-repmat(x2,size(x1,1),1))*(x1-repmat(x2,size(x1,1),1))') / (2*lengthscale^2));
    krr{nx} = @(x) alpha{nx}'*kernel_casadi(X{nx},x,lengthscales(nx,1));
    krr_viz{nx} = @(x) alpha{nx}'*kernel(X{nx},x,lengthscales(nx,1));

end

surrogate_model = @(z) [krr{1}(z); krr{2}(z)];

% you can't because of MATLAB
% you'd have to define 2 krr functions, one for casADi and one for
% visualization...

% % visualizing the learned functions
U = us; 
pts = 225; 
gran1 = (x1_max - x1_min)/(sqrt(pts) - 1);
gran2 = (x2_max - x2_min)/(sqrt(pts) - 1);
[X1, X2] = meshgrid([x1_min:gran1:x1_max], [x2_min:gran2:x2_max]);
Z = [X1(:) X2(:) U*ones(size(X1(:)))];

figure
surf(X1, X2, reshape(krr_viz{1}(Z), size(X1)))
figure
surf(X1, X2, reshape(krr_viz{2}(Z), size(X1)))

%%
% Computing an optimal control sequence

% phase portrait figure
x0 = [1.23; 0.81]; % avg_subopt = 0.0765, avg_opt = 0.0183 (rectangle area per step)
x0 = [1.365; 1.68]; % avg_subopt = 0.0856, avg_opt = 0.0139
x0 = [2.64; 1.33]; % avg_subopt = 0.1391, avg_opt = 0.0201
x0 = [2.7; 0.544]; % avg_subopt = 0.2672, avg_opt = 0.0875

% time-domain simulation
% x0 = [2.7; 0.544];

[x_opti, u_opti] = ocp(surrogate_model, x0, N);

T = [0 N*T_samp];
time = linspace(T(1),T(2),N+1);
x_steps = zeros(N+1,2);
[t,x] = ode45(@(t,x) dynamics(t,x,u_opti), T, x_opti(:,1));
%[t,x] = ode45(@(t,x) dynamics(t,x,us), T, x_opti(:,1));
idx = find(t>=time(1), 1, 'first');
x_true(1,:) = x(idx,:);
for n = 2:N+1
    idx = find(t>=time(n), 1, 'first');
    x_true(n,:) = x(idx,:);
end
figure
plot(x_opti(1,:), x_opti(2,:), '-*');
hold on
plot(x_true(:,1), x_true(:,2), '-*');
plot(xs(1), xs(2), 'kx', 'markersize', 10);
%plot(x(:,1), x(:,2), '-*');
grid on; set(gcf,'color','w');
X_feas = Polyhedron([1 0; -1 0; 0 1; 0 -1],[x_max(1); -x_min(1); x_max(2); -x_min(2)]);
plot(X_feas, 'color', 'red', 'linewidth', 1.5, 'linestyle', '--', 'alpha', 0.02)

%% Part 4
% Collecting data
clear ubs lbs

% extracting a dataset from the system
% D = 600*ones(1,N);
% connected = false; 
% plts = false;
% datasets = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'ocp', plts, surrogate_model);

% TO BE DELETED
D = 200*ones(1,N);
connected = false; 
plts = false;
datasets0 = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'ocp', plts, surrogate_model);
% TO BE DELETED
D = 400*ones(1,N);
connected = false; 
plts = false;
datasets1 = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'ocp', plts, surrogate_model);
% TO BE DELETED
D = 600*ones(1,N);
connected = false; 
plts = false;
datasets2 = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'ocp', plts, surrogate_model);
% TO BE DELETED
D = 800*ones(1,N);
connected = false; 
plts = false;
datasets3 = collect_data(D, @(t,x,u) dynamics(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'ocp', plts, surrogate_model);
% TO BE DELETED

min_queries_dists(datasets, x_opti, u_opti)
min_features_dists(datasets);

%%

% pictures were generated with these
delta_aug = 1.2; 
gamma_aug = 1.2;

% works for the violation example
%delta_aug = 4; 
%gamma_aug = 1.2;

%%

disp([newline 'Computing the sub-optimal bounds:'])

ubs_subopt = zeros(nx,N);
lbs_subopt = zeros(nx,N);
for step = 1:N    
    for state = 1:nx
    
        data = datasets{state,step};
        lengthscale = lengthscales(state,step);
        gamma = gammas(state,step);
        z = [x_opti(:,1)' u_opti(1:step)];
        
        [ub, lb] = approx_bnd(z, data, @(x1,x2) kernel(x1,x2,lengthscale), gamma*gamma_aug, delta_bar*delta_aug);
        ubs_subopt(state,step) = ub;
        lbs_subopt(state,step) = lb;
        
    end
    disp(['Step ' num2str(step) ' done...'])
end

% xx = 1:size(x_opti,2);
% aug_ubs = [x_opti(:,1) ubs];
% aug_lbs = [x_opti(:,1) lbs];

plot_time(x_opti, x_true, ubs_subopt, lbs_subopt, x_min, x_max, 1);
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', gamma augr=' num2str(gamma_aug) ' and aug del=' num2str(delta_aug)])

plot_phase(x_opti, x_true, ubs_subopt, lbs_subopt, x_min, x_max, 2)
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', gamma aug=' num2str(gamma_aug) ' and aug del=' num2str(delta_aug)])
xlim([0.9 3.1]); ylim([0.45 2.05])

plot_errbars(x_opti, x_true, ubs_subopt, lbs_subopt, x_min, x_max, 3);

disp([newline 'Computing the optimal bounds:'])

ubs_opt = zeros(nx,N);
lbs_opt = zeros(nx,N);
for step = 1:N
    for state = 1:nx
    
        data = datasets{state,step};
        lengthscale = lengthscales(state,step);
        gamma = gammas(state,step);
        z = [x_opti(:,1)' u_opti(1:step)];
        
        [ub, lb] = opt_bnd(z, data, @(x1,x2) kernel(x1,x2,lengthscale), gamma*gamma_aug, delta_bar*delta_aug);
        ubs_opt(state,step) = ub;
        lbs_opt(state,step) = lb;
        
    end
    disp(['Step ' num2str(step) ' done...'])
end

% xx = 1:size(x_opti,2);
% aug_ubs = [x_opti(:,1) ubs];
% aug_lbs = [x_opti(:,1) lbs];

plot_time(x_opti, x_true, ubs_opt, lbs_opt, x_min, x_max, 4);
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', gamma augr=' num2str(gamma_aug) ' and aug del=' num2str(delta_aug)])

plot_phase(x_opti, x_true, ubs_opt, lbs_opt, x_min, x_max, 5)
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', gamma augr=' num2str(gamma_aug) ' and aug del=' num2str(delta_aug)])
xlim([0.9 3.1]); ylim([0.45 2.05])

plot_errbars(x_opti, x_true, ubs_opt, lbs_opt, x_min, x_max, 6)


% average area
avg_subopt = sum((ubs_subopt(1,:) - lbs_subopt(1,:)) .* (ubs_subopt(2,:) - lbs_subopt(2,:)))/N;
avg_opt = sum((ubs_opt(1,:) - lbs_opt(1,:)) .* (ubs_opt(2,:) - lbs_opt(2,:)))/N;