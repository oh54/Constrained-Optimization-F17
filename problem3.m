clc; clear all;
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
colorbar

hold on
    fill(x1c1,X2,'b','facealpha',0.2);
    fill(x1c2,X2,'b','facealpha',0.2);
    fill(x1c3,X2,'b','facealpha',0.2);
    axis([0 5 0 3])
hold off

%% sovle the problem

H = [2 0;0 2];
g = [-2;-5];
A = [1 -1 -1 1 0;-2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];

[x,lambda]=EqualityQPSolver(H,g,A,b)

%%
%active set algorithm
% W = {3,5}
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
% W = {5}
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
% W = Ø
A0 = [];
b0 = [];
%solve for W = Ø
x0 = H\-g
%Controle of solution 
c1 = x0(1)-2*x0(2)+2;
c2 = -x0(1)-2*x0(2)+6;
c3 = -x0(1)+2*x0(2)+2;

%%
%active set algorithm
% W = {1}
A1 = [ 1; -2];
b1 = [ -2];
%solve for W = {5}
[x1,lambda1]=EqualityQPSolver(H,g,A1,b1);
%Controle of solution 
c1 = x1(1)-2*x1(2)+2;
c2 = -x1(1)-2*x1(2)+6;
c3 = -x1(1)+2*x1(2)+2;

%%
%Plot the active-set and make a table
figure(1)
hold on
plot(x35(1),x35(2),'*')
plot(x5(1),x5(2),'*')
plot(x0(1),x0(2),'*')
plot(x1(1),x1(2),'*')
hold off
legend('show')
%%
I = eye(size(A,2));
Alin = [A'+I];
f = [1 1 1 1 1];
blin = b';
x=linprog(f,Alin,blin)


