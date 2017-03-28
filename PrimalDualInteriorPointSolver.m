function [x, y, z, s] = PrimalDualInteriorPointSolver(H, g, A, b, C, d)
    [o,p] = size(C);
    [k,l] = size(A);
    [m,n] = size(d);
    
    
    % Initialization of parameters 
    x = [0.5,0.5,0,0,0,0]';
    y = ones(l,1);
    s = ones(m,1);
    z = ones(m,1);
    e = ones(m,1);    
    eta = 0.99;
    max_iter = 100;
    tol = 10^-5;
    
    
    % Compute residuals
    rL = H*x + g - A*y - C*z;
    rA = b - A'*x;
    rC = s + d - C'*x;
    rSZ= z.*s; 
    mu = z'*s/m;

    
    Converged = (norm(rL,inf) <= tol) && ...
                (norm(rA,inf) <= tol) && ...
                (norm(rC,inf) <= tol) && ...
                (abs(mu) <= tol);
    i = 1;

    %%
    while ~Converged && (i<max_iter)
        % LDL factorization of modified KKT        
        H_mod = H + C * diag(z./s) * C';
        KKT = [H_mod, -A; -A', zeros(l,l)];
        [L,D,p] = ldl(KKT,'lower','vector');
        
        %% AFFINE STEP
        % Update residuals
        rL_bar = rL - C * diag(z./s) * (rC - diag(1./z) * rSZ);  

        
        % Calculate affine deltas
        tmp = -[rL_bar;rA];
        X(p) = L'\(D\(L\tmp(p)));
        X = reshape(X, length(X), 1);
        dX_aff = X(1:k);
        dY_aff = X(k+1:k+l);
        dZ_aff = -diag(z./s) * C' * dX_aff + diag(z./s) * (rC - diag(1./z) * rSZ);
        dS_aff = -diag(1./z) * rSZ - diag(s./z) * dZ_aff;
        
        
        % Compute largest alpha_aff
        idxZ = find(dZ_aff < 0.0);
        idxS = find(dS_aff < 0.0);   
        alpha_aff = min([1.0; -z(idxZ,1)./dZ_aff(idxZ,1); -s(idxS,1)./dS_aff(idxS,1)]);

        
        % Dual gap
        mu = z'*s/m;
        mu_aff = (z+alpha_aff*dZ_aff)'*(s+alpha_aff*dS_aff)/m;
        ro = (mu_aff / mu)^3;
        
        
        %% CENTERING STEP
        
        
        % Update residuals
        rSZ_bar= rSZ + dS_aff' * dZ_aff * e - ro * mu * e;
        rL_bar = rL - C * diag(z./s) * (rC - diag(1./z) * rSZ_bar);
        
        
        % Calculate centering deltas
        tmp = -[rL_bar;rA];
        X(p) = L'\(D\(L\tmp(p)));
        X = reshape(X, length(X), 1);
        dX_cen = X(1:k);
        dY_cen = X(k+1:k+l);
        dZ_cen = - diag(z./s) * C' * dX_cen + diag(z./s) *(rC - diag(1./z) * rSZ_bar);
        dS_cen = - diag(1./z) * rSZ_bar - diag(s./z) * dZ_cen;
        
        
        % Compute largest alpha_cen
        idxZ = find(dZ_cen < 0.0);
        idxS = find(dS_cen < 0.0);   
        alpha_cen = min([1.0; -z(idxZ,1)./dZ_cen(idxZ,1); -s(idxS,1)./dS_cen(idxS,1)]);
        alpha = alpha_cen * eta;
        
        
        % Update position
        x = x + alpha * dX_cen;
        y = y + alpha * dY_cen;
        z = z + alpha * dZ_cen;
        s = s + alpha * dS_cen;  
        
        
        % Update residuals
        rL = H*x + g - A*y - C*z;
        rA = b - A'*x;
        rC = s + d - C'*x;
        rSZ= z.*s; 
        mu = z'*s/m;
    
        Converged = (norm(rL,inf) <= tol) && ...
                    (norm(rA,inf) <= tol) && ...
                    (norm(rC,inf) <= tol) && ...
                    (abs(mu) <= tol);
    
        i = i + 1;
    end
    x
    y
    if ~Converged
        x=[];
        y=[];
        z=[];
        s=[];
    end
end