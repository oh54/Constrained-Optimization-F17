function [phi, dphi0,phi0] = phialpha(x,alpha,p,mu)

[f0,g0] = Himmelblau(x);
[c0,~] = constraint_fun(x);
[f,g] = Himmelblau(x+alpha*p);
[c,~] = constraint_fun(x+alpha*p);


phi = f + mu' * abs( min(0, c ));

dphi0 = g0'*p - mu' * abs((min(0,c0)));

phi0 = f0 + mu' * abs((min(0,c0)));
end