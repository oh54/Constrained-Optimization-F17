u_bar = 0.2;
d_0 = 1;
n = 5;

EqualityQPSolverLDL(n,u_bar,d_0)
EqualityQPSolverLU(n,u_bar,d_0)
NullSpaceQR(n,u_bar,d_0)
RangeSpaceSolver(n,u_bar,d_0)