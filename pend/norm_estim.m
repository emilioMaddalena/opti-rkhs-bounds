% Supplementary material for the paper:
% Robust Uncertainty Bounds in Reproducing Kernel Hilbert Spaces:  
% A Convex Optimization Approach'
% Authors: P. Scharnhorst, E. T. Maddalena, Y. Jiang and C. N. Jones
%
% A one-dimensional example

%%
function norm_est = norm_estim(f, kernel, N)
 
    if nargin == 2
        
        N = 5000;
            
    end

    % gathering data
    X = rand(N,1).*(f.xmax-f.xmin) + f.xmin;
    fX = f.foo(X);
    
    % computing the kernel matrix
    jitter = 0.000000001;
    K = kernel(X,X) + jitter*eye(N);
   
    % estimating the RKHS norm from below and applying a safety factor
    norm_est = sqrt((fX'/K)*fX);

end

%EOF