clc; clear all;
X1 = 0:1:10;
X2 = 0:1:10;
[x1,x2] = meshgrid(X1,X2);
%defining the problem
q = (x1-1).^2+(x2-2.5).^2;
x1c1 = 2*X2-2;
x1c2 = -2*X2+6;
x1c3 = 2*X2+2;

%Making the contour plot with the constraints
v = [0:2:10 10:10:100 100:20:200];

contour(x1,x2,q,v, 'linewidth',2)
colorbar

hold on
    fill(X2,x1c1,'b','facealpha',0.2)
    fill(X2,x1c2,'b','facealpha',0.2)
    fill(X2,x1c3,'b','facealpha',0.2)
    axis([0 10 0 10])
hold off

%% sovle the problem

H = [2 0;0 2];
g = [-2;-5];
A = [1 -1 -1 1 0;-2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];

[x,lambda]=EqualityQPSolver(H,g,A,b)

