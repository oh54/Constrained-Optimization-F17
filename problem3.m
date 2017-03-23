clc; clear all; close all;
X1 = 0:1:5;
X2 = 0:1:3;
[x1,x2] = meshgrid(X1,X2);
%defining the problem
q = (x1-1).^2+(x2-2.5).^2;
x1c1 = 2*X2-2;
x1c2 = -2*X2+6;
x1c3 = 2*X2+2;
x1c5 = 0;
x2c5 = 0;

%Making the contour plot with the constraints
v = [0:0.5:3 3:2:15 15:10:100 100:20:200];

figure(1)
contour(x1,x2,q,v, 'linewidth',1);
title('Contour plot of the problem')
xlabel('x_1');
ylabel('x_2');
hold on
    colorbar
    fill(x1c1,X2,'b','facealpha',0.2);
    fill(x1c2,X2,'b','facealpha',0.2);
    fill(x1c3,X2,'b','facealpha',0.2);
    axis([0 5 0 3]);
hold off

%% sovle the problem

H = [2 0;0 2];
g = [-2;-5];
A = [1 -1 -1 1 0;-2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];

[x,lambda]=EqualityQPSolver(H,g,A,b)

%%
%active set algorithm
% first workingset W = {3,5}
A35 = [-1  0;2 1];
b35 = [-2; 0];
%solve for W = {3,5}
[x35,lambda35]=EqualityQPSolver(H,g,A35,b35);
%Controle of solution 
c1 = x35(1)-2*x35(2)+2;
c2 = -x35(1)-2*x35(2)+6;
c3 = -x35(1)+2*x35(2)+2;

%%
%active set algorithm
%second worikngset W = {5}
A5 = [ 0; 1];
b5 = [ 0];
%solve for W = {5}
[x5,lambda5]=EqualityQPSolver(H,g,A5,b5);
%Controle of solution 
c1 = x5(1)-2*x5(2)+2;
c2 = -x5(1)-2*x5(2)+6;
c3 = -x5(1)+2*x5(2)+2;

%%
%active set algorithm
% thrid workingset W = Ø
A0 = [];
b0 = [];
%solve for W = Ø
x0 = H\-g;
lambda0 = 0;
%Controle of solution 
c1 = x0(1)-2*x0(2)+2;
c2 = -x0(1)-2*x0(2)+6;
c3 = -x0(1)+2*x0(2)+2;

%%
%active set algorithm
% fouth workingset W = {1}
A1 = [ 1; -2];
b1 = [ -2];
%solve for W = {1}
[x1,lambda1]=EqualityQPSolver(H,g,A1,b1);
%Controle of solution 
c1 = x1(1)-2*x1(2)+2;
c2 = -x1(1)-2*x1(2)+6;
c3 = -x1(1)+2*x1(2)+2;


%%
%Plot the active-set and make a table
figure(1)
hold on
p1 = plot(x35(1),x35(2),'r*', 'MarkerSize',15);
p2 = plot(x5(1),x5(2),'b*', 'MarkerSize',15);
p3 = plot(x0(1),x0(2),'g*', 'MarkerSize',15);
p4 = plot(x1(1),x1(2),'k*', 'MarkerSize',15);
hold off
legend([p1 p2 p3 p4],'Workingset 1','Workingset 2','Workingset 3', 'Workingset 4');
%%
Name = {'Working set 1';'Working set 2';'Working set 3';'Working set 4'};
Active_cons = [35; 5; 2; 1];
x_1 = [x35(1); x5(1); x0(1); x1(1)];
x_2 = [x35(2); x5(2); x0(2); x1(2)];
lambda_1 = [lambda35(1); lambda5; NaN; lambda1];
lambda_2 = [lambda35(2); NaN; NaN; NaN];

T = table(Active_cons,x_1,x_2,lambda_1, lambda_2,'RowNames',Name)
%%
I = eye(size(A,2));
Alin = [A' I];
f = [0 0 1 1 1 1 1]';
blin = b';
x_lp=linprog(f,-Alin,-blin,[],[],zeros(7,1),inf)
figure(1)
hold on
    p5 = plot(x_lp(1),x_lp(2),'m*','MarkerSize',15);
hold off
legend([p5],'linprog point')
%% 
x0 = x_lp;
x_qp = quadprog(H,g,-A',-blin,[],[],zeros(2,1),inf,x0)


