function [ x, lambda ] = EqualityQPSolver(H,g,A,b)
    % K matrix
    
    sizes = size(A);
    ncols = sizes(2);
    
    KKT_A = [ H -A ; -A' zeros(ncols, ncols) ];
    KKT_b = -[ g ; b ];

    [L,D,p] = ldl(KKT_A, 'lower', 'vector');
    KKT_x = L'\( D \(L\KKT_b(p)));
    
    x = KKT_x(1:end-2);
    lambda = KKT_x(end-1:end);
    
end

