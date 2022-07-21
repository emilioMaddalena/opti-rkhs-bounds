function [x,u] = gen_feats(N, D, xmin, xmax, umin, umax, connected)
% Generate features (random x and u) in a safe way, that is, with a minimum
% distance among datapoints

    disp([newline 'Generating features:'])

    threshold = 1e-1;

    nx = numel(xmin);
    nu = numel(umin);
    
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

