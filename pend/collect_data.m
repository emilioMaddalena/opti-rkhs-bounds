function DATASET = collect_data(N, D, dynamics, xmin, xmax, umin, umax, delta_bar, connected, plts)
%
% Returns a {x_dim, N} cell array.
%
% Each cell is the data matrix associated with a system state dimension and
% a certain time-step. In other words, each datum (row of the matrix)
% contains the current state, the control sequence, and a (noise corrupted)
% component of the next state.
%
% The first cell is a [D, x_dim+(u_dim)+1] matrix.
% The second cell is a [D, x_dim+(u_dim*2)+1] matrix.
% ...
%  
    global T_samp
    
    plot_aug = 1.5;
    
    if ~exist('connected')
        connected = true;
    end
    
    if ~exist('plts')
        plts = false;
    end

    if connected
    % One datum at a time, simulate for the longest horizon and 
    % place each segment into its dataset
    
        if size(D,2) ~= 1, error('Connected data requires D to be scalar...'); end
        
        % OLD: x0 = (xmax - xmin).*rand(2,D) + xmin;
        % OLD: u  = (umax - umin).*rand(1,N,D) + umin;
        
        [x0,u] = gen_feats(N, D, xmin, xmax, umin, umax, connected);
        
        T = [0 N*T_samp];
        time = linspace(T(1),T(2),N+1);
        x_steps = zeros(N+1,2);

        for d = 1:D

            % Integrating the system
            [t,x] = ode45(@(t,x) dynamics(t,x,u(1,:,d)), T, x0(:,d));

            % Extracting x0
            idx = find(t>=time(1), 1, 'first');
            x_steps(1,:) = x(idx,:);

            for n = 2:N+1

                idx = find(t>=time(n), 1, 'first');
                x_steps(n,:) = x(idx,:);

                if plts 
                    figure(n-1)
                    plot(x_steps(1,1),x_steps(1,2),'ko','markersize',8,'linewidth',2); hold on; grid on
                    plot(x_steps(n,1),x_steps(n,2),'kx','markersize',8,'linewidth',2);
                    plot(x_steps(1:n,1), x_steps(1:n,2),'-b','linewidth',2);
                    axis(plot_aug*[xmin(1) xmax(1) xmin(2) xmax(2)])
                    title(['Data collection for ' num2str(n-1) ' control moves!'])
                end
                
                % adding noise
                del = rand(2,1) * 2*delta_bar - delta_bar;

                DATASET{1,n-1}(d,:) = [x_steps(1,:) u(1,1:n-1,d) (x_steps(n,1) + del(1))];
                DATASET{2,n-1}(d,:) = [x_steps(1,:) u(1,1:n-1,d) (x_steps(n,2) + del(2))];
            end
        end
        
     elseif ~connected
     % Fixes the horizon, build entire dataset for it, increases horizon
     
        if size(D,2) == 1, D = D*ones(N,1); end
                
        for n = 1:N
    
            T = [0 n*T_samp];
            time = linspace(T(1),T(2),n+1);
            x_steps = zeros(n,2);

            [X0,U] = gen_feats(n, D(n), xmin, xmax, umin, umax, connected);
            
            for d = 1:D(n)
                
                % OLD: x0 = (xmax - xmin).*rand(2,1) + xmin;
                % OLD: u  = (umax - umin).*rand(1,n) + umin;
                
                x0 = X0(:,d);
                u  = U(d,1:n);

                [t,x] = ode45(@(t,x) dynamics(t,x,u), T, x0);

                for i = 1:n+1
                    idx = find(t>=time(i), 1, 'first');
                    x_steps(i,:) = x(idx,:);
                end

                if plts
                    figure(n)
                    plot(x_steps(1,1),x_steps(1,2),'ko','markersize',5,'linewidth',1.5); hold on; grid on 
                    plot(x_steps(end,1),x_steps(end,2),'kx','markersize',5,'linewidth',1.5);
                    plot(x_steps(:,1), x_steps(:,2),'-b','linewidth',1)
                    axis(plot_aug *[xmin(1) xmax(1) xmin(2) xmax(2)])
                    title(['Data collection for ' num2str(n) ' control moves!'])
                end
                
                % adding noise
                del = rand(2,1) * 2*delta_bar - delta_bar;
                
                DATASET{1,n}(d,:) = [x_steps(1,:) u (x_steps(end,1) + del(1))];
                DATASET{2,n}(d,:) = [x_steps(1,:) u (x_steps(end,2) + del(2))];
            end

        end
    end
     
    disp('Done collecting data!')
end