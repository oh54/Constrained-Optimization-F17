% problem 4 - Markowitz
H = [ 2.30 0.93 0.62 0.74 -0.23 ; 0.93 1.40 0.22 0.56 0.26 ; 0.62 0.22 1.80 0.78 -0.27 ; 0.74 0.56 0.78 3.40 -0.56 ;-0.23 0.26 -0.27 -0.56 2.60 ];
returns = [ 15.10 ; 12.50 ; 14.70 ; 9.02 ; 17.68 ];

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
ub = [];
% Q2
% assuming not shortsales maximum return 17.68 and minimal 9.02

% Q3
[x,fval,exitflag,output,lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub);

% optimal portfolio weights
x
% corresponding risk 
fval


% Q4
% do the above computation many times and plot risk-return