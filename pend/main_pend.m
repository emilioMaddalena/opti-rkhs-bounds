% Supplementary material for the paper:
% "
%
%

clear all
clc

addpath('./utils')

rng(1)
global T_samp 

constants

kernel = @(x1,x2,lengthscale) exp(-dist(x1,x2').^2 / (2*lengthscale^2));

% PART1

% extracting a dataset from the system
N = 3;  
D = [500 500 500];
delta_bar = 0; 
connected = false; 
visuals = false;
dataset = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'grid', connected, visuals);

% estimate the kernel lengthscales
lengthscale_range = 1:0.2:3;
folds = 10;
lambda = 0.01;
[lengthscales, losses] = cv(folds, kernel, lambda, dataset, lengthscale_range);

% estimate ground-truths RKHS norms
gammas = estimate_rkhs(dataset, @(x1,x2,lengthscale) kernel(x1,x2,lengthscale), lengthscales);

%%
% PART2
% Learning a 1-step ahead KRR model for the system

% extracting a dataset from the system
N = 1;  
D = 100; 
delta_bar = 0.01; 
connected = false; 
visuals = false;    
dataset = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'grid', connected, visuals);

for nx = 1:nx
    
    X{nx} = dataset{nx,1}(:,1:end-1);
    y{nx} = dataset{nx,1}(:,end);

    jitter = 1e-8;
    n_data = size(X{nx},1);
    K{nx} = kernel(X{nx},X{nx},lengthscales(nx,1)) + jitter*eye(n_data);
    alpha{nx} = (K{nx} + n_data*lambda*eye(n_data))\y{nx}; 
    
    kernel_casadi = @(x1,x2,lengthscale) exp(-diag((x1-repmat(x2,size(x1,1),1))*(x1-repmat(x2,size(x1,1),1))') / (2*lengthscale^2));
    krr{nx} = @(x) alpha{nx}'*kernel_casadi(X{nx},x,lengthscales(nx,1));

end

surrogate_model = @(z) [krr{1}(z); krr{2}(z)];

% you can't because of MATLAB
% you'd have to define 2 krr functions, one for casADi and one for
% visualization...

% % visualizing the learned functions
% U = 0; 
% P = 225; 
% gran1 = (x1_max - x1_min)/(sqrt(P) - 1);
% gran2 = (x2_max - x2_min)/(sqrt(P) - 1);
% [X1, X2] = meshgrid([x1_min:gran1:x1_max], [x2_min:gran2:x2_max]);
% Z = [X1(:) X2(:) U*ones(size(X1(:)))];
% 
% surf(X1, X2, reshape(krr{1}(Z), size(X1)))
% figure
% surf(X1, X2, reshape(krr{2}(Z), size(X1)))

% Computing an optimal control sequence

x0 = [pi/2; 0];
N = 3;

[x_opti, u_opti] = ocp(surrogate_model, x0, N);

% plot(x_opti(1,:),x_opti(2,:))
% figure 
% plot(u_opti)

%% Part 4
% Building an optimal certification box

clear ubs lbs

% extracting a dataset from the system
N = 3;  
D = [50 700 1050]; % 1200
delta_bar = 0.01; 
connected = false; 
visuals = false;

rng(1)
datasets = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, 'grid', connected, visuals);

min_queries_dists(datasets, x_opti, u_opti)
min_features_dists(datasets);

delta_bar = 0.03;
aug_factor = 1.4;

disp([newline 'Computing the optimal bounds:'])
for step = 1:N
    
    for state = 1:nx
    
        data = datasets{state,step};
        lengthscale = lengthscales(state,step);
        gamma = gammas(state,step)*aug_factor;
        z = [x_opti(:,1)' u_opti(1:step)];
        
        [ub, lb] = opt_bnd(z, data, @(x1,x2) kernel(x1,x2,lengthscale), gamma, delta_bar);
        ubs(state,step) = ub;
        lbs(state,step) = lb;
        
    end

    disp(['Step ' num2str(step) ' done...'])

end

xx = 1:size(x_opti,2);
aug_ubs = [x_opti(:,1) ubs];
aug_lbs = [x_opti(:,1) lbs];

figure
for state = 1:nx

    subplot(1,2,state)
    plot(xx, x_opti(state,:), 'k-', 'marker', 'o', 'linewidth', 2); hold on
    fill([xx fliplr(xx)], [aug_lbs(state,:) fliplr(aug_ubs(state,:))], 'y', 'facealpha', 0.15, 'linewidth', 1.5);
    
    plot(xx, x_max(state)*ones(size(xx)), 'k--', 'linewidth', 2)
    plot(xx, x_min(state)*ones(size(xx)), 'k--', 'linewidth', 2)
    %ylim([x_min(state)*1.5 x_max(state)*1.5]);
    xlabel('time step'); xticks(1:numel(xx));
    ylabel(['x' num2str(state)]);
    
end
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', aug factor=' num2str(aug_factor) ' and del bar=' num2str(delta_bar)])

ubs