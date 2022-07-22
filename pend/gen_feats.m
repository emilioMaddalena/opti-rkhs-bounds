function [x,u] = gen_feats(N, D, xmin, xmax, umin, umax, method, connected)
% Generate features (random x and u) in a safe way, that is, with a minimum
% distance among datapoints

    threshold = 1e-1;

    nx = numel(xmin);
    nu = numel(umin);
    
    if strcmp(method,'rnd')

        d = 1;
        z = zeros(D,nx+(nu*N));

        while d <= D

            x_candidate = (xmax - xmin).*rand(nx,1) + xmin;
            u_candidate = (umax - umin).*rand(nu,N,1) + umin;

            z_candidate = [x_candidate' u_candidate];

            if min(dist(z,z_candidate')) >= threshold
                z(d,:) = z_candidate;
                d = d + 1;
                if mod(d,100) == 0, disp([num2str(d) ' samples collected...']); end
            end

        end
        
    elseif strcmp(method,'grid')
        
        %D = 500;
        d_dim = ceil(D^(1/(nx+N*nu))) - 1;

        x1_gran = (xmax(1) - xmin(1)) / d_dim;
        x2_gran = (xmax(2) - xmin(2)) / d_dim;
        u_gran = (umax - umin) / d_dim;

        x1_vec = xmin(1):x1_gran:xmax(1);
        x2_vec = xmin(2):x2_gran:xmax(2);
        u_vec = umin:u_gran:umax;
        
        % replicate the control vec N times
        coords = {x1_vec, x2_vec};
        for i = 1:N
            coords{nx+i} = u_vec;
        end
        
        %[X1,X2,U] = ndgrid(coords{:});
        %X1 = X1(:); 
        %X2 = X2(:);
        %U = U(:);
        %z = [X1(:) X2(:) U(:)];
        
        data = cell(1,nx+N*nu);
        [data{:}] = ndgrid(coords{:});
        z = [];
        for i = 1:(nx+N*nu)
            z = [z data{i}(:)];
        end   
        
        % dropping rnd elements to return exactly D elements
        for i = 1:(size(z,1)-D)
          
            idx = randi([1,size(z,1)]);
            z(idx,:) = [];
            
        end
        
    end

    % Final formatting
    if connected

        % weird formatting...
        x = z(:,1:nx)';
        for i = 1:D
            u(:,:,i) = z(i,nx+1:end);
        end

    else

        x = z(:,1:nx)';
        u = z(:,nx+1:end);

    end

    
end

