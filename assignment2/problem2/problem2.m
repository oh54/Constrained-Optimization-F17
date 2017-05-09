clear;

x_0 = [ -1.8 1.7 1.9 -0.8 -0.8 ]';

lambda_0 = [ 1 1 1 ]';
x_star = [ -1.71 1.59 1.82 -0.763 -0.763 ]';

x_true = [-1.7171 1.5957 1.8272 -0.7636 -0.7636 ]';

maxiter = 10;
epsilon = 10^-4;

[x1, fval1, lambda1, k1, X1] = EqSQP(x_0, lambda_0, maxiter, epsilon);

[x2, fval2, lambda2, k2, X2] = EqSQP_BFGS(x_0, lambda_0, maxiter, epsilon);

% parameter for scaling alpha in linesearch
tau = 0.9;
% parameter for scaling sufficient decrease term in linesearch
eta = 0.4;
[x3, fval3, lambda3, k3, X3, alphas] = EqSQP_BFGS_linesearch(x_0, lambda_0, maxiter, epsilon, tau, eta);

% hold on
% X = X1;
% [m,n] = size(X);
% e_n = X(:,1:n-1) - repmat(x_true, 1, n-1);
% e_nplus1 = X(:,2:n) - repmat(x_true, 1, n-1);
% e_n = sqrt(sum(e_n.^2, 1));
% e_nplus1 = sqrt(sum(e_nplus1.^2, 1));
% plot(-log(e_n), -log(e_nplus1), 'Color', 'b');
% xlabel('-log(e_k)');
% ylabel('-log(e_{k+1})');
X = X2;
[m,n] = size(X);
e_n = X(:,1:n-1) - repmat(x_true, 1, n-1);
e_nplus1 = X(:,2:n) - repmat(x_true, 1, n-1);
e_n = sqrt(sum(e_n.^2, 1));
e_nplus1 = sqrt(sum(e_nplus1.^2, 1));
p1 = plot(-log(e_n), -log(e_nplus1), '--mo', 'Color', 'b', 'LineWidth', 2);
hold on
X = X3;
[m,n] = size(X);
e_n = X(:,1:n-1) - repmat(x_true, 1, n-1);
e_nplus1 = X(:,2:n) - repmat(x_true, 1, n-1);
e_n = sqrt(sum(e_n.^2, 1));
e_nplus1 = sqrt(sum(e_nplus1.^2, 1));
p2 = plot(-log(e_n), -log(e_nplus1), '-.mo', 'Color', 'r', 'LineWidth', 2, 'MarkerSize',3);
xlabel('-log(e_k)');
ylabel('-log(e_{k+1})');
legend([p1 p2],{'Algo2','Algo3'});
hold off


%d1 = digits(3);
%latex(vpa(sym(X3')))
%digits(d1);

