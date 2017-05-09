function [x,lambda]=EqualityQPSolverLU(H,g,A,b)
KKT_A = [H -A; A' zeros(length(b),length(b))];
KKT_b = [-g;b];
[L,U,p] = lu(KKT_A,'vector');
KKT_x = U \(L\ KKT_b(p) );
x = KKT_x(1:length(H));
lambda = KKT_x((length(H)+1):length(KKT_x));
end


