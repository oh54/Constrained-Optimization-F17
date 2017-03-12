function [x, lambda] = EqualityQPSolverLDL(n, u_bar, d_0)
    [KKT_A, KKT_b, H, g, A, b] = ConstructKKT(n, u_bar, d_0);
    n = size(A,1);
    m = size(A,2);
    
    [L,D,p] = ldl(KKT_A,'lower','vector');
    KKT_x(p) = L'\(D\(L\KKT_b(p)));
    KKT_x = KKT_x';
    
    x = KKT_x(1:n);
    lambda = KKT_x(n+1:n+m); 
end
