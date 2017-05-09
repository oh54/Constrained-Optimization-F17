function [ x_k, fval, lambda_k, k, X, delta ] = IneqSQP_BFGS_trust( x_0, lambda_0, maxiter, epsilon, gamma, eta )
% Inequality constrained non-linear programming solver using
% sequential quadratic programming approach with damped BFGS approximation
% to the hessian of the lagrangian and trust region
x_k = x_0;
X = [ x_k ];
lambda_k = lambda_0;
k = 0;
B_k = eye(2);
grad_fk = gradient_f(x_k);
A_k = jacobian_c(x_k);
delta_k = 0.5;
delta = [ delta_k ];
mu_k = 1;

while (k < maxiter) 
    k
    c_k = c(x_k);
    %delta_k
    % KKT conditions satisfied
    if max(abs(grad_fk - A_k' * lambda_k)) < epsilon && max(abs(c_k)) >= -epsilon && ...
            min(lambda_k) >= -epsilon && max(abs(c_k .* lambda_k)) < epsilon
        break
    end
    H_quad = [ B_k zeros(2) ; zeros(2) zeros(2) ];
    f_quad = [grad_fk ; mu_k ; mu_k ];
    A_quad = [-A_k eye(2) ; 1 0 0 0 ; -1 0 0 0 ; 0 1 0 0 ; 0 -1 0 0];
    b_quad = [c_k ; delta_k ; delta_k ; delta_k ; delta_k];
    lb = [ -Inf -Inf 0 0 ]';
    [p_t, ~, ~, ~, lambda_t] = quadprog(H_quad, f_quad, A_quad, b_quad, [], [], lb, []);
    p_k = p_t(1:2);
    lambda_hat = lambda_t.ineqlin(1:2);
   
    % p.550
    ared_k = phi_1(x_k, mu_k) - phi_1(x_k + p_k, mu_k);
    q_zero = q_mu(f(x_k), grad_fk, B_k, mu_k, [ 0 0 ]', c_k, A_k);
    q_p = q_mu(f(x_k), grad_fk, B_k, mu_k, p_k, c_k, A_k);
    pred_k = q_zero - q_p;
    rho_k = ared_k / pred_k;
    
    p_lambda = lambda_hat - lambda_k;
    lambda_kplus1 = lambda_k + p_lambda;
    
    if rho_k > eta
        x_kplus1 = x_k + p_k;
        delta_kplus1 = delta_k;
    else
        x_kplus1 = x_k;
        delta_kplus1 = gamma * norm(p_k, inf);
        
        x_k = x_kplus1;
        delta_k = delta_kplus1;
        X = [X x_k];
        delta = [delta delta_k];
        lambda_k = lambda_kplus1;
        k = k + 1;  
        continue
    end
    
    % p536 18.13
    s_k = x_kplus1 - x_k;
    grad_fkplus1 = gradient_f(x_kplus1);
    A_kplus1 = jacobian_c(x_kplus1);
    grad_Lkplus1 = grad_fkplus1 - A_kplus1' * lambda_kplus1;
    grad_Lk = grad_fk - A_k'*lambda_kplus1;
    y_k = grad_Lkplus1 - grad_Lk;
    
    % p537 
    theta_k = get_theta_k(B_k, s_k, y_k);
    r_k = theta_k * y_k + (1-theta_k) * B_k * s_k;
     
    % standard BFGS update with y_k replaced by r_k
    % standard BFGS update would not necessarily retain the positive
    % semidefiniteness of the B_k matrix
    add1 = (B_k * s_k * s_k' * B_k) / (s_k' * B_k * s_k);
    add2 = (r_k * r_k') / (s_k' * r_k);
    B_kplus1 = B_k - add1 + add2;
    
    % update everything for next iteration
    B_k = B_kplus1;
    grad_fk = grad_fkplus1;
    A_k = A_kplus1;
    x_k = x_kplus1;
    delta_k = delta_kplus1;
    X = [X x_k];
    delta = [delta delta_k];
    lambda_k = lambda_kplus1;
    k = k + 1;  
end
fval = f(x_k);