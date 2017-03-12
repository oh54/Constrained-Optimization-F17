u_bar = 0.2;
d_0 = 1;
n = 5;
clc
[x, lambda] = EqualityQPSolverLDL(n,u_bar,d_0)
[x, lambda] = EqualityQPSolverLU(n,u_bar,d_0)
[x, lambda] = NullSpaceQR(n,u_bar,d_0)
[x, lambda] = RangeSpaceSolver(n,u_bar,d_0)