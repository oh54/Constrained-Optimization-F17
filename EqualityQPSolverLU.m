function [x, lambda] = EqualityQPSolverLU(n, u_bar, d_0)
    [KKT_A, KKT_b, H, g, A, b] = ConstructKKT(n, u_bar, d_0);
    n = size(A,1);
    m = size(A,2);
    
    [L,U,p]=lu(KKT_A,'vector');
    
    KKT_x = U\(L\KKT_b(p));
    x = KKT_x(1:n);
    lambda = KKT_x(n+1:n+m);
end
