clc; clear all; close all;

x_0 = [ 4 3 ]';
%x_0 = [ -4 0 ]';
%x_0 = [ 1,-1 ]';

lambda_0 = [ 1 1 ]';

%maxiter = 20;
maxiter = 15;
epsilon = 10^-1;

%[x_bfgs, fval_bfgs, lambda_bfgs, k1, X_bfgs] = IneqSQP_BFGS(x_0, lambda_0, maxiter, epsilon);

tau = 0.9;
eta = 0.4;

%[x_line, fval_line, lambda_line, k2, X_line, alphas] = IneqSQP_BFGS_linesearch(x_0, lambda_0, maxiter, epsilon, tau, eta);

gamma = 0.5;
[x_trust, fval_trust, lambda_trust, k3, X_trust, delta_trust] = IneqSQP_BFGS_trust(x_0, lambda_0, maxiter, epsilon, gamma, eta);

X = X_trust;

d1 = digits(3);
latex(vpa(sym(X')))
digits(d1);


step = 0.05;
x1vec = -5:step:5;
x2vec = -5:step:5;
[x1,x2] = meshgrid(x1vec, x2vec);

obj = (x1.^2 + x2 - 11).^2 + (x1 + x2.^2 -7).^2;

c1 = (x1 + 2).^2 - x2 >= 0;
c2 = -4*x1 + 10*x2 >= 0;

% % Q1 
v = [0:2:6 6:6:24 24:20:124 124:24:244];
contour(x1,x2,obj,v, 'linewidth',1);
xlabel('x_1');
ylabel('x_2');
colorbar
hold on
calpha=0.8;
msize=0.15;
scatter(x1(~c1), x2(~c1), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
scatter(x1(~c2), x2(~c2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
axis([-5 5 -5 5]);
%hold off

msize = 10;
calpha = 1;

i = 1
[m, n] = size(X);
while i <= n
    hold on;
    scatter(X(1,i), X(2,i), msize, 'red', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
    if i < 13
        text(X(1,i) + 0.1, X(2,i), int2str(i-1));  
    end
    
    rectangle('Position', [X(1,i)-delta_trust(i) X(2,i)-delta_trust(i) delta_trust(i)*2 delta_trust(i)*2]);
    
    i = i + 1;
end

