% problem 4 - Markowitz
cov = [ 2.30 0.93 0.62 0.74 -0.23 ; 0.93 1.40 0.22 0.56 0.26 ; 0.62 0.22 1.80 0.78 -0.27 ; 0.74 0.56 0.78 3.40 -0.56 ;-0.23 0.26 -0.27 -0.56 2.60 ];
returns = [ 15.10 ; 12.50 ; 14.70 ; 9.02 ; 17.68 ];

H = cov .* 2;
%H = cov;
% Q1
% objective function should minimize risk
% there are two equality constraints, one so that return = 10 and the other
% so the proportions of money invested sum to 1
% covariance matrix is essentially the Hessian matrix and x is the
% vector with the weights in each asset
% if shortsales are not allowed, there's additional inequality constraint
% that all weights must be larger/equal to zero

f = [ 0 ; 0 ; 0 ; 0 ; 0 ];
A = -returns';
b = -10;
Aeq = [ 1 1 1 1 1 ];
beq = 1;
lb = [ 0 ; 0 ; 0 ; 0 ; 0 ];
ub = [ 1 ; 1 ; 1 ; 1 ; 1 ];
% Q2
% assuming not shortsales maximum return 17.68 and minimal 9.02

% Q3
%[x,fval,exitflag,output,lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
% R = 10 portfolio
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_10 = x;
fval_10 = fval;

% R = 17
b = -17;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_17 = x;
fval_17 = fval;

b = -17.5;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_175 = x;
fval_175 = fval;


% min variance portfolio
b = 0;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);
x_min = x;
fval_min = fval;

% Q4
plot_variances = [];
plot_returns = [];
for R = 9:0.05:17.70
    b = -R;
    [x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);    

    plot_returns = [plot_returns, x' * returns];
    plot_variances = [plot_variances, fval];

end
p1 = plot(plot_variances, plot_returns, 'LineWidth', 2);
ylim([8, 21]);
xlim([0, 4]);
title('Efficient frontier');
xlabel('Risk (portfolio variance)');
ylabel('Return (%)');
hold on;
%scatter(fval_10, x_10' * returns, 20); 

%txt1 = '\leftarrow Minimum variance portfolio';
%text(fval_10 + 0.05, x_10' * returns + 0.03, txt1)

% Q5
returns = [ 15.10 ; 12.50 ; 14.70 ; 9.02 ; 17.68 ; 2.0 ];
% add a column and row of zeros to the covariance matrix, representing that the riskfree
% security isn't correlated with any of the previous securities
cov = [cov, zeros(5,1)];
cov = [cov ; zeros(1,6)];
H = 2 .* cov;

f = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0];
lb = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0];
ub = [ 1 ; 1 ; 1 ; 1 ; 1 ; 1];

% Q6
Aeq = [ 1 1 1 1 1 1 ];
beq = 1;
A = -returns';

plot_variances = [];
%plot_returns = zeros(21, 1);
plot_returns = [];
for R = 7:0.05:17.70
    b = -R;
    [x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);    

    plot_returns = [plot_returns, x' * returns];
    plot_variances = [plot_variances, fval];

end

hold on;
p2 = plot(plot_variances, plot_returns, 'LineWidth', 2);

hold on;
vars = diag(cov);

for i = 1:5
   s = scatter(vars(i), returns(i), 20); 
   txt = int2str(i);
   text(vars(i)+0.1, returns(i), txt);
    
end


legend([p1 p2],{'Initial','With risk free asset'});

% Q7
b = -15;
[x,fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub);  
x_15_rf = x;
fval_15_rf = fval;
r_15_rf = returns' * x_15_rf;

hold on;
p_15 = scatter(fval_15_rf, r_15_rf, 20); 
txt15 = '\leftarrow P_{15}';
text(fval_15_rf + 0.1, r_15_rf -0.05, txt15);


% Problem 5 question 4
H = H;
g = f;
C = [-A' eye(6)];
A = Aeq';
b = beq;
d = [-15 ; zeros(6,1)];

