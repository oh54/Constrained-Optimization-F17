clear;

x_0 = [ -1.8 1.7 1.9 -0.8 -0.8 ]';
x_star = [ -1.71 1.59 1.82 -0.763 -0.763 ]';
lambda_0 = [ 10 10 10 ];

x_k = x_0;
lambda_k = lambda_0;

k = 0;

while (k < 10)
    f_k = f(x_k);
    grad_f = gradient_f(x_k);

    hessian_L = hessian_f(x_k) - (lambda_k(1) * hessian_ci(x_k,1) + ...
                                  lambda_k(2) * hessian_ci(x_k,2) + ...
                                  lambda_k(3) * hessian_ci(x_k,3));
    c_k = c(x_k);
    A_k = jacobian_c(x_k);
    
    % quadprog(hessian_L, grad_f, [], [], A_k, -c_k)
    [p, lambda] = EQPSolver(hessian_L, grad_f, A_k', -c_k);
    
    x_k = x_k + p;
    lambda_k = lambda;
    k = k + 1;
end

