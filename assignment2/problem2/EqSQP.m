function [ x_k, fval, lambda_k, k, X ] = EqSQP( x_0, lambda_0, maxiter, epsilon )
% Equality constrained non-linear programming solver using
% sequential quadratic programming approach


x_k = x_0;
X = [x_k];
lambda_k = lambda_0;
k = 0;

while (k < maxiter)
    grad_fk = gradient_f(x_k);

    hessian_L = hessian_f(x_k) - (lambda_k(1) * hessian_ci(x_k,1) + ...
                                  lambda_k(2) * hessian_ci(x_k,2) + ...
                                  lambda_k(3) * hessian_ci(x_k,3));
    c_k = c(x_k);
    A_k = jacobian_c(x_k);
    
    % 8 element vector (n+m)
    F_k = [ grad_fk - A_k' * lambda_k ; c_k ];
    % first order KKT conditions 
    if max(abs(F_k)) < epsilon
        break
    end
    
    % quadprog(hessian_L, grad_fk, [], [], A_k, -c_k)
    [p, lambda] = EQPSolver(hessian_L, grad_fk, A_k', -c_k);
    
    x_k = x_k + p;
    X = [X x_k];
    lambda_k = lambda;
    k = k + 1;
    
end
fval = f(x_k);

end

