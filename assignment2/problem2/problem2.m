clear;

x_0 = [ -1.8 1.7 1.9 -0.8 -0.8 ]';
lambda_0 = [ 10 10 10 ]';
x_star = [ -1.71 1.59 1.82 -0.763 -0.763 ]';

maxiter = 20;
epsilon = 10^-3;

[x1, fval1, lambda1, k1, X1] = EqSQP(x_0, lambda_0, maxiter, epsilon);

[x2, fval2, lambda2, k2, X2] = EqSQP_BFGS(x_0, lambda_0, maxiter, epsilon);

% parameter for scaling alpha in linesearch
tau = 0.9;
% parameter for scaling sufficient decrease term in linesearch
eta = 0.1;
[x3, fval3, lambda3, k3, X3] = EqSQP_BFGS_linesearch(x_0, lambda_0, maxiter, epsilon, tau, eta);


% x_k = x_0;
% lambda_k = lambda_0;
% instead of calculating hessian exactly, some approximation could be used
%hessian_L = hessian_f(x_k) - (lambda_k(1) * hessian_ci(x_k,1) + ...
%                              lambda_k(2) * hessian_ci(x_k,2) + ...
%                              lambda_k(3) * hessian_ci(x_k,3));
%B_k = hessian_L;    
% B_k = eye(5);
% grad_fk = gradient_f(x_k);
% A_k = jacobian_c(x_k);
% c_k = c(x_k);
% 
% [p_k, lambda_hat] = EQPSolver(B_k, grad_fk, A_k', -c_k);
% p_lambda = lambda_hat - lambda_k;
% sigma = 1;
% rho = 0.5;
% mu_k = 0;
% margin = 0.01;
% mu_k = get_mu(mu_k, grad_fk, p_k, B_k, c_k, sigma, rho, margin);
% 
% eta = 0.25;
% alpha_k = 1;
% 
% tau = 0.5;
% tau_alpha = 0.25;
% 
% i = 0;
% while phi_1(x_k + alpha_k * p_k, mu_k) > phi_1(x_k, mu_k) + get_decrease_term(mu_k, eta, alpha_k, grad_fk, p_k, c_k) 
%     alpha_k = tau_alpha * alpha_k;
%     
%     i = i+1;
% end


