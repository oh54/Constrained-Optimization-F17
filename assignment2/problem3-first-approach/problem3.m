
clc; clear all; close all;
%difining the problem 
step = 0.1;
x1vec = -5:step:5;
x2vec = -5:step:5;
[x1,x2] = meshgrid(x1vec, x2vec);

q = (x1.^2+x2-11).^2+(x1+x2.^2-7).^2;
c1 = (x1+2).^2-x2 >= 0;
c2 = -4*x1+10*x2 >= 0;

%defining constant etc. 
n = 2;
x0 = [4;3];
B0 = eye(n);
x = x0;
B_k = B0;
tol = 10^-4;
i = 0;
theta = 1; 
alpha = 1;
mu = [0;0];
x_true = [3;2];

%making the contour plot of the problem
v = [0:2:6,6:6:24,24:20:124,124:24:244];
figure(1); hold on;
    title('contour plot of the problem');
    xlabel('x_1');
    ylabel('x_2');
    colorbar
contour(x1,x2,q,v);
    calpha=0.5;
    msize=0.6;
    scatter(x1(~c1), x2(~c1), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
    scatter(x1(~c2), x2(~c2), msize, 'black', 'MarkerFaceAlpha', calpha, 'MarkerEdgeAlpha', calpha);
    axis([-5 5 -5 5]);

scatter(x0(1),x0(2),[10],'r','filled');
    text(x0(1),x0(2),'S','color','r','FontSize',10,'HorizontalAlignment','right');

for k = 1:1;
    [~,nablaf_k] = Himmelblau(x);
    g = nablaf_k;
        if norm(g,inf) < tol
            break 
        else            
            H = B_k;
            [ck,Ak] = constraint_fun(x);
        end
    [p,~,~,~,lambda] = quadprog(H,g,-Ak',ck);
   
    s = lambda.ineqlin;
    mu = max(abs(s),0.5*(mu+abs(s)));
    %alpha = linesearch(x, alpha, p, mu) %the line search part
    p = alpha*p; % if the algorithm should be used without line search algorithm, alpha can be redefine to equal 1 below the lineseach function, or the line with the linesearch comand can be marked as a comment.    
    [gradc] = Himmelblau_L_grad(x);
    c = gradc;
    nablax_L = g - c*s;
    
    x = x + p;
    
    [~,nablaf_k_new] = Himmelblau(x);
    [gradc] = Himmelblau_L_grad(x);
    c = gradc;
    y = nablaf_k_new - c*s - nablax_L;
 
       
    if p'*y >= 1/5*p'*B_k*p; 
        theta = 1;
    else
        theta = 4/5* (p'*B_k*p)/(p'*B_k*p-p'*y);
    end
    r_k = theta*y+(1-theta)*B_k*p;
    
    B_k = B_k - (B_k*p*p'*B_k)/(p'*B_k*p)+(r_k*r_k')/(r_k'*p);
    i = i +1;
    
    scatter(x(1),x(2),10,'black','filled');
    text(x(1),x(2),num2str(i),'FontSize',9);
    
    num = [i (x(1)) (x(2)) (((x(1)-x_true(1))^2+(x(2)-x_true(2))^2)^(1/2)) (nablax_L(1)) (nablax_L(2))];
    T(i,:)=num;
end
x_end = x;
B_end = B_k;
scatter(x_end(1),x_end(2),[15],'R','filled');
text(x_end(1),x_end(2),'F','color','r','FontSize',10,'HorizontalAlignment','right');
hold off
Table = array2table(T,'VariableNames',{'iterationNumber' 'x1' 'x2' 'error' 'Nablax_L_x1' 'Nablax_L_x2'});


%%
%the trust region 





