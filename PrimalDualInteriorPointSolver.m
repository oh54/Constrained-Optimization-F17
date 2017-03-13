function [x, lambda] = PrimalDualInteriorPointSolver(H, g, A, b, C, d, x, y, Z, S)
    m = size(s,1);
    e = ones(m,1);
    k = size(A,1);
    l = size(A,2);
    eta = 0.995;
    
    max_iter = 1000;
    eps = 10^-5;
    
    % Compute residuals
    rL = H*x + g - A*y - C*Z;
    rA = b - A'*x;
    rC = s + d - C'*x;
    rSZ= Z*S*e; 
    mu = z'*s/m;
            
    i = 1;
    
    % Get norms of residuals
    norm_L = norm(rL,inf);
    norm_A = norm(rA,inf);
    norm_C = norm(rC,inf);
    
    while i<max_iter & norm_L <= eps*max(1,normL) & norm_C <= eps*max(1,norm_C) & norm_A <= eps*max(1,norm_A) & abs(mu) > eps_mu   
              
        % LDL factorization of modified KKT
        inv_S = pinv(S);
        inv_Z = pinv(Z);
        H = H + C*(inv_S*Z)*C';
        KKT = [H, -A; -A', zeros(l,l)];
        [L,D,p] = ldl(KKT,'lower','vector');

        
        %% AFFINE STEP
        % Update residuals
        rL_bar = rL - C * (inv_S*Z) * (rC-inv_Z*rSZ);
        
        % Calculate affine deltas
        tmp = -[rL_bar;rA];
        X(p) = L'\(D\(L\tmp(p)));
        dX_aff = X(1:k);
        dY_aff = X(k+1:k+l);
        
        tmp = inv_S*Z
        dZ_aff = -tmp * C' * dX_aff + tmp * (rC - inv_Z * rSZ);
        dS_aff = -inv_Z * rSC - inv_Z * S * dZ_aff;
        
        % Compute largest alpha_aff
        alpha_aff = 1;
        % >> To be implemented <<

        % Dual gap
        mu = z'*s/m;
        mu_aff = (z+alpha_aff*dZ_aff)'*(s+alpha_aff*dS_aff)/m;
        ro = (mu_aff/mu)^3;
        
        
        %% CENTERING STEP
        
        % Update residuals
        rSZ_bar= rSZ + dS_aff * dZ_aff*e - ro*mu*e;
        rL_bar = rL - C * (inv_S*Z) * (rC-inv_Z*rSZ_bar);
        
        % Calculate centering deltas
        tmp = -[rL_bar;rA];
        X(p) = L'\(D\(L\tmp(p)));
        dX_cen = X(1:k);
        dY_cen = X(k+1:k+l);
        dZ_cen = - (inv_S*Z) * C' * dX_cen + (inv_S*Z)*(rC - inv_Z * rSZ_bar);
        dS_cen = - inv_Z * rSC_bar - inv_Z * S * dZ_cen;
        
        
        % Compute largest alpha_cen
        alpha_cen = 1;
        % >> To be implemented <<
        
        alpha = alpha_cen * eta;
               
        % Update position
        x = x + alpha * dX_cen;
        y = y + alpha * dY_cen;
        z = z + alpha * dZ_cen;
        s = s + alpha * dS_cen;
                
        % Update residuals
        rL = H*x + g - A*y - C*Z;
        rA = b - A'*x;
        rC = s + d - C'*x;
        rSZ= Z*S*e; 
        
        % Get norms of residuals
        norm_L = norm(rL,inf);
        norm_A = norm(rA,inf);
        norm_C = norm(rC,inf);
    
        i = i + 1;
    end
    
    x = 0;
    lambda = 0;
end