function [x, lambda] = PrimalDualInteriorPointSolver(H, g, A, b, C, d, x, y, Z, S)
    m = size(s,1);
    e = ones(m,1);
    k = size(A,1);
    l = size(A,2);

    while 1
              
        % LDL factorization of modified KKT
        inv_S = pinv(S);
        inv_Z = pinv(Z);
        H = H + C*(inv_S*Z)*C';
        KKT = [H, -A; -A', zeros(l,l)];
        [L,D,p] = ldl(KKT,'lower','vector');

        % Compute residuals
        rL = H*x + g - A*y - C*z;
        rA = b - A'*x;
        rC = s + d - C'*x;
        rSZ= Z*S*e; 

        % Calculate deltas
        tmp = -[rL;rA];
        X(p) = L'\(D\(L\tmp(p)));
        dX_aff = X(1:k);
        dY_aff = X(k+1:k+l);
        
        tmp = inv_S*Z
        dZ_aff = -tmp * C' * dX_aff + tmp * (rC - inv_Z * rSZ);
        dS_aff = -inv_Z * rSC - inv_Z * S * dZ_aff;
        
        % Compute largest alpha_affine
      
        
        
        
        
    end
    
    x = 0;
    lambda = 0;
end