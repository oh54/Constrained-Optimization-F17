function [x,lambda]=EQPSolver(H,g,A,b)
F = [H -A; A' zeros(length(b),length(b))];
d = [-g;b];
%LU factorisation
[L,U,p] = lu(F,'vector');
%solve F*z = d
z = U \(L\ d(p) );
%separate x and lambda
x = z(1:length(H));
lambda = z((length(H)+1):length(z));
end

