function [x, lambda] = RangeSpaceSolver(n, u_bar, d_0)
    % This solver is already sparse as the Cholesky factorization with
    % 3 ouput values require sparse input
    
    [H,g,A,b] = ConstructFactorization(n,u_bar,d_0);
    
    n = size(A,1);
    m = size(A,2);  
    
    % Hv = g
    [L,~,s] = chol(sparse(H),'lower','vector');
    v = L'\(L\g(s));
    
    % inv(H) ~ H
    H_A = A'*H*A;
    
    % H_A*lambda = b + A'*v
    [L,~,s] = chol(sparse(H_A),'lower','vector');
    tmp = b + A'*v;
    lambda(s) = L'\(L\tmp(s));
    lambda = lambda';

    % H*x = A*lambda - g
    tmp = A*lambda - g;
    x = H \ tmp;
end