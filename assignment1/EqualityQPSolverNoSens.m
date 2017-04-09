function [x, lambda, dx, dlambda] = EqualityQPSolverNoSens(H,g,A,b)
    n = size(A,1);
    m = size(A,2);  
    KKT_A = [H, -A; -A', zeros(m,m)];
    KKT_b = -[g;b];

    [L,D,p] = ldl(KKT_A,'lower','vector');
    KKT_x(p) = L'\(D\(L\KKT_b(p)));
    KKT_x = KKT_x';
    %the final resulat 
    x = KKT_x(1:n);
    lambda = KKT_x(n+1:n+m);
    %the sensitivity solver
    %out_sens = -KKT_A^(-1);
    %dx = out_sens(:,1:length(g));
    %dlambda = out_sens(:,length(g)+1:end);   
end