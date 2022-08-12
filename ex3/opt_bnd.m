function [ub,lb] = opt_bnd(x, data, kernel, gamma, del_bar)
%OPT_BND Summary of this function goes here
%   Detailed explanation goes here

    function bound = compute_bound(x, dir)

        X = data(:,1:end-1);
        y = data(:,end);
        N = size(X,1);
        
        if N <= 300, jitter = 1e-5; else, jitter = 1e-5; end
        K = kernel(X,X) + jitter*eye(N); 

        c  = sdpvar(N,1);
        cz = sdpvar(1);

        if strcmp(dir,'up'), obj = -cz; elseif strcmp(dir,'dw'), obj = cz; end
        const = [[c; cz]'*([K kernel(X,x); kernel(X,x)' kernel(x,x)]\[c; cz]) <= gamma^2];
        const = [const, norm(c-y, inf) <= del_bar];

        options = sdpsettings('verbose', 0, 'solver', 'gurobi');
        sol = optimize(const, obj, options);

        if sol.problem ~= 0
            disp(sol.problem)
            error('Opt problem is infeasible!');
        end
        
        bound = value(cz);

    end

    ub = compute_bound(x, 'up');
    lb = compute_bound(x, 'dw');

end
