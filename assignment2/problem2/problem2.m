clear;

x_0 = [ -1.8 1.7 1.9 -0.8 -0.8 ]';
lambda_0 = [ 10 10 10 ]';
x_star = [ -1.71 1.59 1.82 -0.763 -0.763 ]';

maxiter = 10;
epsilon = 10^-9;

[x1, fval1, lambda1, k1] = EqSQP(x_0, lambda_0, maxiter, epsilon);

[x2, fval2, lambda2, k2] = EqSQP_BFGS(x_0, lambda_0, maxiter, epsilon);

