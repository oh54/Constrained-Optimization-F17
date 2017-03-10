function [KKT_A, KKT_b, H, g, A, b] = ConstructKKT(n, u_bar, d_0)
    [H,g,A,b] = ConstructFactorization(n,u_bar,d_0);
    
    n = size(A,1);
    m = size(A,2);  
    KKT_A = [H, -A; -A', zeros(m,m)];
    KKT_b = -[g;b];
end