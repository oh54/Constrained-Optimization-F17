function alpha = linesearch(x, alpha, p, mu)
c1= 0.1;
tau = 0.5;
[phi,dphi0,phi0] = phialpha(x,alpha,p,mu);
iter = 1;
max_iter = 100;
while (phi > phi0 + c1*alpha*dphi0) && (iter<=max_iter)%The Armijo condition
    
    alpha = alpha*tau;
    [phi, dphi0,phi0] = phialpha(x,alpha,p,mu);
    iter = iter + 1;
end
alpha = alpha;
end