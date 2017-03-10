function [x, lambda] = RangeSpaceSolver(n, u_bar, d_0)
    [H,g,A,b] = ConstructFactorization(n,u_bar,d_0);
    
    n = size(A,1);
    m = size(A,2);  
    
    L = chol(H,'lower');
    v = H\g;
    
    ro = 1/(v'*g);
    H_A = ( eye(n,n) - ro * v * g') * ( eye(n,n) - ro * g * v' ) + ro * v * v'; 
    L_A = chol(H_A,'lower');
    
    
    lambda = H_A \ ( b + A'*v);
    x = H \ (A*lambda - g);
end