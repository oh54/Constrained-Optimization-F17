% problem 4 - Markowitz
H = [ 2.30 0.93 0.62 0.74 -0.23 ; 0.93 1.40 0.22 0.56 0.26 ; 0.62 0.22 1.80 0.78 -0.27 ; 0.74 0.56 0.78 3.40 -0.56 ;-0.23 0.26 -0.27 -0.56 2.60 ];
returns = [ 15.10 ; 12.50 ; 14.70 ; 9.02 ; 17.68 ];

H = H .* 2;

% Q1
% objective function should minimize risk
% there are two equality constraints, one so that return = 10 and the other
% so the proportions of money invested sum to 1
% covariance matrix is essentially the Hessian matrix and x is the
% vector with the weights in each asset
% if shortsales are not allowed, there's additional inequality constraint
% that all weights must be larger/equal to zero

f = [ 0 ; 0 ; 0 ; 0 ; 0 ];
A = [];
b = [];
Aeq = [ 1 1 1 1 1 ; returns' ];
beq = [ 1 ; 10 ];
lb = [ 0 ; 0 ; 0 ; 0 ; 0 ];
ub = [ 1 ; 1 ; 1 ; 1 ; 1 ];
% Q2
% assuming not shortsales maximum return 17.68 and minimal 9.02

% Q3
%[x,fval,exitflag,output,lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
% R = 10 portfolio
beq = [ 1 ; 17 ];
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_17 = x;
fval_17 = fval;

% R = 17
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_10 = x;
fval_10 = fval;


% min variance portfolio
Aeq = [ 1 1 1 1 1 ];
beq = 1;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_min = x;
fval_min = fval;

% Q4
% The efficient frontier can be traced out as a combinations of any two
% efficient portfolios

%plot_variances = zeros(21, 1);
plot_variances = [];
%plot_returns = zeros(21, 1);
plot_returns = [];
%i = 0
for alpha = -2:0.05:2
    %i = i+1;
    beta = 1 - alpha;
    xvec = x_17 .* alpha + x_min .* beta;
    
    %plot_returns(i) = returns' * xvec;
    %plot_variances(i) = xvec' * H * xvec;
    
    plot_returns = [plot_returns, returns' * xvec];
    plot_variances = [plot_variances, xvec' * H * xvec];
end

figure();
plot(plot_variances, plot_returns);
xlim([0, 8]);
hold on;

% Q5
H = [ 2.30 0.93 0.62 0.74 -0.23 ; 0.93 1.40 0.22 0.56 0.26 ; 0.62 0.22 1.80 0.78 -0.27 ; 0.74 0.56 0.78 3.40 -0.56 ;-0.23 0.26 -0.27 -0.56 2.60 ];
returns = [ 15.10 ; 12.50 ; 14.70 ; 9.02 ; 17.68 ; 2.0 ];

% add a column and row of zeros to H, representing that the riskfree
% security isn't correlated with any of the previous securities
H = [H, zeros(5,1)];
H = [H ; zeros(1,6)];

f = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0];
lb = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0];
ub = [ 1 ; 1 ; 1 ; 1 ; 1 ; 1];

% Q6

Aeq = [ 1 1 1 1 1 1 ];
beq = 1;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_min_rf = x;
fval_min_rf = fval;


Aeq = [ 1 1 1 1 1 1 ; returns' ];
beq = [ 1 ; 10 ];
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_10_rf = x;
fval_10_rf = fval;


plot_variances = [];
%plot_returns = zeros(21, 1);
plot_returns = [];
%i = 0
for alpha = -100:0.05:100
    %i = i+1;
    beta = 1 - alpha;
    xvec = x_10_rf .* alpha + x_min_rf .* beta;
    
    %plot_returns(i) = returns' * xvec;
    %plot_variances(i) = xvec' * H * xvec;
    
    plot_returns = [plot_returns, returns' * xvec];
    plot_variances = [plot_variances, xvec' * H * xvec];
end

plot(plot_variances, plot_returns, 'g');
%xlim([0, 8]);

% Q7
beq = [ 1 ; 15 ];
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_15_rf = x;
fval_15_rf = fval;


