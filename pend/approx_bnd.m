function [ub,lb] = approx_bnd(x, data, kernel, gamma, del_bar)
%Description here

    function bound = compute_bound(x, dir)
        
        % build nominal KRR model
        jitter = 1e-5;
        lambda = 1e-2;

        X = data(:,1:end-1);
        y = data(:,end);

        n_data = size(X,1);
        K = kernel(X,X) + jitter*eye(n_data);
        alpha = (K + n_data*lambda*eye(n_data))\y; 
        krr = @(x) alpha'*kernel(X,x);

%         % compute Delta
%         delta = sdpvar(n_data,1);
%         constr = -del_bar <= delta <= del_bar;
%         cost = (delta'/K)*delta - 2*(y'/K)*delta;
%         temp = optimize(constr, cost, sdpsettings('verbose', 0, 'solver', 'gurobi'));
%         delta = value(delta);
%         Delta = -(delta'/K)*delta + 2*(y'/K)*delta;
% 
%         % compute the indivdual terms and the bounds
%         s = @(x) y'*(K\kernel(X,x));
%         s_norm = sqrt(y'*(K\y));
%         P = @(x) sqrt(real(diag(kernel(x,x) - (kernel(x,X)/K)*kernel(X,x))));
%         p = @(x) (del_bar * vecnorm((kernel(x,X)/K)', 1))';
%         %q = @(x) (y' * ((K+(1/(N*lam))*K*K)\kernel(X,x)))';
%         q = @(x) abs(s(x) - krr(x));
% 
%         S = @(x) P(x).*sqrt(gamma^2 + Delta - s_norm^2) + abs(p(x)) + q(x);

        % compute Delta
        nu = sdpvar(n_data,1);
        cost = 1/4*(nu'*K*nu) + nu'*y + del_bar*norm(nu,1);
        temp = optimize([], cost, sdpsettings('verbose', 0, 'solver', 'gurobi'));
        nu = value(nu);
        Delta = value(cost);

        % compute the indivdual terms and the bounds
        s = @(x) y'*(K\kernel(X,x));
        s_norm = sqrt(y'*(K\y));
        P = @(x) sqrt(real(diag(kernel(x,x) - (kernel(x,X)/K)*kernel(X,x))));
        p = @(x) (del_bar * vecnorm((kernel(x,X)/K)', 1))';
        q = @(x) abs(s(x) - krr(x));

        S = @(x) P(x).*sqrt(gamma^2 + Delta) + abs(p(x)) + q(x);

        if strcmp(dir,'up'), bound = krr(x) + S(x);
        elseif strcmp(dir,'dw'), bound = krr(x) - S(x); end
    
    end
    
    ub = compute_bound(x, 'up');
    lb = compute_bound(x, 'dw');
    
end