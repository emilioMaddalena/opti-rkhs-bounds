function [gammas, lengthscales] = estimate_rkhs(dataset, kernel, lengthscale_range)
%
% Returns two (x_dim, N) arrays: 
% gammas - the estimated RKHS norms for each x_dim and time step
% lengthscales - the estimated best hyperparameters
%
% TODO: apply safety factor

    jitter = 1e-5;

    nx = size(dataset, 1); % number of states
    nu = size(dataset{1,1}, 2) - (nx+1); % number of control inputs
    nsteps = size(dataset, 2);  % num of steps 
    
    for time_step = 1:nsteps

        for x_dim = 1:nx
        
            data = dataset{x_dim,time_step};
            norm_ests = [];
            
            for lengthscale = lengthscale_range 

                % gathering data
                N = size(data,1);
                X = data(:,1:end-1);
                y = data(:,end);

                % computing the kernel matrix
                K = kernel(X,X,lengthscale) + jitter*eye(N);

                % estimating the RKHS norm from below 
                norm_est = sqrt((y'/K)*y);
                norm_ests = [norm_ests norm_est];

            end
            
            lowest_norm = min(norm_ests);
            
%             p = plot(lengthscale_range, norm_ests, 'linewidth', 2); grid on; xlabel('lengthscale'); ylabel('RKHS estimate');
%             disp(['Lowest norm: ', num2str(lowest_norm)])
%             pause
%             close

            % save lowest RKHS norm and its associated lengthscale
            gammas(x_dim,time_step) = lowest_norm;
            lengthscales(x_dim,time_step) = lengthscale_range(norm_ests == lowest_norm);
            
        end
        
    end

end

