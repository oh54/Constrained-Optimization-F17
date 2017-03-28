clc; clear all; close all;
step = 0.05;
x1vec = -1:step:5;
x2vec = -1:step:5;
[x1,x2] = meshgrid(x1vec, x2vec);

q = (x1-1).^2 + (x2-2.5).^2;

c1 = x1 - 2*x2 + 2 >= 0;
c2 = -x1 - 2*x2 + 6 >= 0;
c3 = -x1 + 2*x2 + 2 >= 0;
c4 = x1 >= 0;
c5 = x2 >= 0;

% % Q1 
% v = [0:0.5:3 3:2:15 15:10:100 100:20:200];
% contour(x1,x2,q,v, 'linewidth',1);
% title('Feasible starting point obtained with linprog')
% xlabel('x_1');
% ylabel('x_2');
% colorbar
% hold on
% calpha=0.6;
% msize=0.05;
% scatter(x1(~c1), x2(~c1), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c2), x2(~c2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c3), x2(~c3), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c4), x2(~c4), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c5), x2(~c5), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% axis([-1 5 -1 5]);
% hold off

% Q2 - KKT conditions
H = [2 0;0 2];
g = [-2;-5];
A = [1 -1 -1 1 0;-2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];

% % Q3
% 
% 
% % Q4
%[x,lambda]=EqualityQPSolver(H,g,A,b);


% Q5 - active set
% W_1 = {3,5};
x_0 = [2 ; 0];
A35 = [-1  0 ; 2 1];
b35 = [-2; 0];
[x1,lambda1]=EqualityQPSolver(H,g,A35,b35);
% testing substituting for lambda
lambda1_2 = A35 \ (H * x_0 + g);
p1 = x1 - x_0;
x_1 = x_0 + p1;


% dropped constraint 3 because lambda_3 = -2
% W_2 = {5}
A5 = [ 0; 1];
b5 = 0;
[x2,lambda2]=EqualityQPSolver(H,g,A5,b5);
p2 = x2 - x_1;
x_2 = x_1 + p2;


% dropped constraint 5 because lambda_5 = -5
% W_3 = {}
A0 = [];
b0 = [];
x3 = H\-g;
p3 = x3 - x_2;
lambda3 = 0;
% TODO change to calculation
p3 = 0.6 * p3;
x_3 = x_2 + p3;


% W_4 = {1}
A1 = [ 1; -2];
b1 = -2;
[x4,lambda4]=EqualityQPSolver(H,g,A1,b1);
p4 = x4 - x_3;
x_4 = x_3 + p4;

% W_5 = {1}
A1 = [ 1; -2];
b1 = -2;
[x5,lambda5]=EqualityQPSolver(H,g,A1,b1);
p5 = x5 - x_4;
x_5 = x_4 + p5;

% initial feasible point
% x_0
% % x AFTER each iterations
% x_1
% x_2
% x_3
% x_4
% x_5
% 
% % steps on each iteration
% p1
% p2
% p3
% p4
% p5
%
% % lambdas on each iteration
% lambda1
% lambda2
% lambda3
% lambda4
% lambda5
% 
% msize = 20;
% calpha = 1;
% hold on;
% scatter(x_1(1),x_1(2), msize, 'red', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% text(x_1(1) + 0.1,x_1(2), '0 & 1');
% scatter(x_2(1),x_2(2), msize, 'red', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% text(x_2(1) + 0.1,x_2(2), '2');
% scatter(x_3(1),x_3(2), msize, 'red', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% text(x_3(1) + 0.1,x_3(2), '3');
% scatter(x_4(1),x_4(2), msize, 'red', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% text(x_4(1) + 0.1,x_4(2), '4 & 5');

% Q6 - comment on Langrange multipliers
% in report

% Q7 - explain active set method
% in report

% Q8
f = [];
[x_lp,fval_lp,exitflag_lp,output_lp,lambda_lp] = linprog(f, -A',-b);
% feasible starting point
x_lp;
% all lambdas zero (up to accuracy), no active constraints
lambda_lp.ineqlin;

% msize = 10;
% calpha = 1;
% hold on;
% scatter(x_lp(1),x_lp(2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% text(x_lp(1) + 0.1,x_lp(2) + 0.02, '\leftarrow Feasible starting point');

% Q9

% specify starting conditions
H = [2 0;0 2];
g = [-2;-5];
A = [1 -1 -1 1 0;-2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];
maxiter = 10;
% linear program starting point
%W = [];
%x_k = x_lp;
% starting point in book
W = [3 5];
x_k = [2;0];
% some infeasible starting point 
%W = [];
%x_k = [1; 4];

% big M starting conditions
% A_M = [A ; ones(1,5)];
% A_M = [A_M [0 ; 0 ; 1]];
% b_M = [b ; 0];
% H_M = zeros(3);
% H_M(1,1) = 2;
% H_M(2,2) = 2;
% g_M = [g ; 9999];
% 
% H = H_M;
% g = g_M;
% A = A_M;
% b = b_M;
% maxiter = 10;
% % book starting point
% W = [3,5,6];
% x_k = [2; 0 ; 0];
% % % some infeasible starting point
% % W = [6];
% %x_k = [1; 4 ; 0];


for k = 1:maxiter
    k
    %x_k
    %W
    if isempty(W)
        xpos = H\-g;
        lambda = [];
    else
        [xpos,lambda]=EqualityQPSolverNoSens(H,g,A(:,W),b(W));
    end
    % unscaled step p_k
    lambda
    x_k
    p_k = xpos - x_k;
    p_k
    
    if all(p_k == 0)
        if (all(lambda >= 0))
            % zero step and lambdas non-negative, solution found
            break 
        else
            % zero step but exists negative lambda, remove it and cont.
            [minlambda, j] = min(lambda);
            W = W(W ~= W(j));
        end
    else
        % Compute alpha
        % indexes not in W
        notW = setdiff(1:5, W);
        % logical vector of indexes not in W that break feasibility
        logBnotW = logical(A(:,notW)' * xpos < b(notW));
        % indexes not in W that break feasibility
        breakers = notW(logBnotW);
        
        % p.469 16.41
        numerator = b(breakers) - A(:,breakers)' * x_k;
        denominator = A(:,breakers)' * p_k;
        alpha_k =  min(1, min(numerator ./ denominator));
        
        if isempty(alpha_k)
            alpha_k = 1;
        end
        alpha_k
        x_k = x_k + alpha_k * p_k;
        
        % if there are blocking constraints
        if length(breakers) >= 1
            % add one of them to W
            W = sort([W breakers(1)]);
        end
    end
end

% Q9
% big M method for quadprog
A_M = [A' ones(5,1)];
H_M = zeros(3);
H_M(1,1) = 2;
H_M(2,2) = 2;
g_M = [g ; 9999];
lb = [-inf -inf 0];
[x_Mq,fval_Mq,exitflag_Mq,output_Mq,lambda_Mq] = quadprog(H_M, g_M, -A_M, -b, [], [], lb, []);

x_Mq
