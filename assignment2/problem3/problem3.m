clc; clear all; close all;
%difining the problem
%making a contour plot of the problem 
step = 0.1;
x1vec = -5:step:5;
x2vec = -5:step:5;
[x1,x2] = meshgrid(x1vec, x2vec);

q = (x1.^2+x2-11).^2+(x1+x2.^2-7).^2;
c1 = (x1+2).^2-x2 >= 0;
c2 = -4*x1+10*x2 >= 0;

v = [0:2:6,6:6:24,24:20:124,124:24:244];
contour(x1,x2,q,v);
title('contour plot of the problem')
xlabel('x_1');
ylabel('x_2');
colorbar
hold on
 calpha=0.5;
 msize=0.6;
 scatter(x1(~c1), x2(~c1), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
 scatter(x1(~c2), x2(~c2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
 axis([-5 5 -5 5]);
hold off
%%
n = 2;
x0 = [2;2];
B0 = eye(n);
x = x0;
B_k = B0;
tol = 10^-4;
for k = 1:1000
    [~,nablaf_k] = Himmelblau(x);
    g = nablaf_k;
        if norm(g,inf) < tol
            break 
        else            
            H = B_k;
            [ck,Ak] = constraint_fun(x);
    [p,~,~,~,lambda] = quadprog(H,g,-Ak,ck);
    s = lambda.ineqlin;
    x = x + p;
    [~,nablaf_k_new] = Himmelblau(x);
    y = nablaf_k_new - nablaf_k;
    
    B_k = B_k - (B_k*s*s'*B_k)/(s'*B_k*s)+(y*y')/(y'*s);
end
x_end = x
B_end = B_k
