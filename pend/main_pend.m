clear all
clc

rng(1)
global T_samp 
constants

kernel = @(x1,x2,lengthscale) exp(-dist(x1,x2').^2 / (2*lengthscale^2));
 
% D = 500;
% d_dim = ceil(D^(1/3));
% 
% x1_gran = (x1_max - x1_min) / d_dim;
% x2_gran = (x2_max - x2_min) / d_dim;
% u_gran = (u_max - u_min) / d_dim;
% 
% x1_vec = x1_min:x1_gran:x1_max;
% x2_vec = x2_min:x2_gran:x2_max;
% u_vec = u_min:u_gran:u_max;
% 
% [X1,X2,U] = ndgrid(x1_vec,x1_vec,u_vec);


% PART1
% Estimating the system's RKHS norms 

% extracting a dataset from the system
N = 3;  
D = [500 500 500];
delta_bar = 0; 
connected = false; 
visuals = false;
  
dataset = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, connected, visuals);

% % estimating RKHS norms and kernel lengthscales
% lengthscale_range = 3:0.1:5;
% [gammas, lengthscales] = estimate_rkhs(dataset, @(x1,x2,lengthscale) kernel(x1,x2,lengthscale), lengthscale_range);

lengthscales = [3 1 1; 0.8 2 0.9];
gammas = estimate_rkhs2(dataset, @(x1,x2,lengthscale) kernel(x1,x2,lengthscale), lengthscales);

%% PART2
% Learning a 1-step ahead KRR model for the system

% extracting a dataset from the system
N = 1;  
D = 100; 
delta_bar = 0.01; 
connected = false; 
visuals = false;    

dataset = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, connected, visuals);

for nx = 1:nx
    
    X{nx} = dataset{nx,1}(:,1:end-1);
    y{nx} = dataset{nx,1}(:,end);

    lam = 1e-4;
    jitter = 1e-8;
    n_data = size(X{nx},1);
    K{nx} = kernel(X{nx},X{nx},lengthscales(nx,1)) + jitter*eye(n_data);
    alpha{nx} = (K{nx} + n_data*lam*eye(n_data))\y{nx}; 
    
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

x0 = [1; 0];
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
D = [200 400 500]; 
delta_bar = 0.01; 
connected = false; 
visuals = false;

dataset = collect_data(N, D, @(t,x,u) pend(t,x,u), x_min, x_max, u_min, u_max, delta_bar, connected, visuals);

% Study the distances among dataset features
min_dist = inf(1,N);
for step = 1:N
    data = dataset{1,step}(:,1:end-1);
    d = D(step);
    for i = 1:d
        data_temp = data;
        sample = data_temp(i,:);
        data_temp(i,:) = [];
        sample_min_dist = min(dist(data_temp,sample'));
        if sample_min_dist < min_dist(step), min_dist(step) = sample_min_dist; end
    end
end
min_dist

% Study the distance of the dataset and the query points
for step = 1:N
    disp(['Step ' num2str(step) ':'])
    z = [x_opti(:,1)' u_opti(1:step)];
    Z = dataset{1,step}(:,1:end-1);
    disp(['Min distance = ' num2str(min(dist(Z,z')))]);
end


disp([newline 'Computing the optimal bounds:'])
aug_factor = 1;
for step = 1:N
    
    for state = 1:nx
    
        data = dataset{state,step};
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
    ylim([x_min(state)*1.5 x_max(state)*1.5])
    xlabel('time step'); xticks(1:numel(xx));
    ylabel(['x' num2str(state)]);
    
end
sgtitle(['Nominal traj and optimal bounds for D=' num2str(D) ', aug factor=' num2str(aug_factor) ' and del bar=' num2str(delta_bar)])
