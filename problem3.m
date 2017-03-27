clc; clear all; close all;
step = 0.05;
x1vec = -4:step:6;
x2vec = -2.5:step:7.5;
[x1,x2] = meshgrid(x1vec, x2vec);

q = (x1-1).^2 + (x2-2.5).^2;

c1 = x1 - 2*x2 + 2 >= 0;
c2 = -x1 - 2*x2 + 6 >= 0;
c3 = -x1 + 2*x2 + 2 >= 0;
c4 = x1 >= 0;
c5 = x2 >= 0;

% Q1 
% v = [0:0.5:3 3:2:15 15:10:100 100:20:200];
% contour(x1,x2,q,v, 'linewidth',1);
% title('Contour plot of the problem')
% xlabel('x_1');
% ylabel('x_2');
% colorbar
% hold on
% calpha=0.6;
% msize=0.01;
% scatter(x1(~c1), x2(~c1), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c2), x2(~c2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c3), x2(~c3), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c4), x2(~c4), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% scatter(x1(~c5), x2(~c5), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
% axis([-4 6 -2.5 7.5]);
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
x_1 = [2 ; 0];
A35 = [-1  0 ; 2 1];
b35 = [-2; 0];
[x35,lambda35]=EqualityQPSolver(H,g,A35,b35);
% testing substituting for lambda
lambda35_2 = A35 \ (H * x_1 + g);
x_2 = x_1;

% dropped constraint 3 because lambda_3 = -2
% W_2 = {5}
A5 = [ 0; 1];
b5 = 0;
[x5,lambda5]=EqualityQPSolver(H,g,A5,b5);
p2 = x5 - x_2;
x_3 = x_2 + p2;

% dropped constraint 5 because lambda_5 = -5
% W_3 = {}
A0 = [];
b0 = [];
xe = H\-g;
p3 = xe - x_3;
lambdae = 0;
% TODO change to calculation
alpha3 = 0.6;
x_4 = x_3 + alpha3 * p3;


% W_4 = {1}
A1 = [ 1; -2];
b1 = -2;
[x1,lambda1]=EqualityQPSolver(H,g,A1,b1);
p4 = x1 - x_4;
x_5 = x_4 + p4;


% W_5 = {1}
A1 = [ 1; -2];
b1 = -2;
[x1,lambda1]=EqualityQPSolver(H,g,A1,b1);
p5 = x1 - x_5;
x_5 = x_5 + p4;



% %%
% %Plot the active-set and make a table
% figure(1)
% hold on
% p1 = plot(x35(1),x35(2),'r*', 'MarkerSize',15);
% p2 = plot(x5(1),x5(2),'b*', 'MarkerSize',15);
% p3 = plot(x0(1),x0(2),'g*', 'MarkerSize',15);
% p4 = plot(x1(1),x1(2),'k*', 'MarkerSize',15);
% hold off
% legend([p1 p2 p3 p4],'Workingset 1','Workingset 2','Workingset 3', 'Workingset 4');
% %%
% Name = {'Working set 1';'Working set 2';'Working set 3';'Working set 4'};
% Active_cons = [35; 5; 2; 1];
% x_1 = [x35(1); x5(1); x0(1); x1(1)];
% x_2 = [x35(2); x5(2); x0(2); x1(2)];
% lambda_1 = [lambda35(1); lambda5; NaN; lambda1];
% lambda_2 = [lambda35(2); NaN; NaN; NaN];
% 
% T = table(Active_cons,x_1,x_2,lambda_1, lambda_2,'RowNames',Name)
% %%
% I = eye(size(A,2));
% Alin = [A' I];
% f = [0 0 1 1 1 1 1]';
% blin = b';
% x_lp=linprog(f,-Alin,-blin,[],[],zeros(7,1),inf)
% figure(1)
% hold on
%     p5 = plot(x_lp(1),x_lp(2),'m*','MarkerSize',15);
% hold off
% legend([p5],'linprog point')
% %% 
% x0 = x_lp;
% x_qp = quadprog(H,g,-A',-blin,[],[],zeros(2,1),inf,x0)
% 
% 
