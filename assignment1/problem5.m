u_bar = 0.2;
d_0 = 1;
n = 5;
clc
[x, lambda] = EqualityQPSolverLDL(n,u_bar,d_0)

[KKT_A, KKT_b, H, g, A, b] = ConstructKKT(n, u_bar, d_0);
[m,n] = size(A);
C = zeros(m,n);
d = zeros(n,1);

[x, y, z, s] = PrimalDualInteriorPointSolver(H, g, A, b, C, d)


