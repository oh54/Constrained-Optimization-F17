function [x, y, z, s] = PrimalDualInteriorPointSolver(H, g, A, b, C, d)
    [o,p] = size(C);
    [k,l] = size(A);
    [q,r] = size(d)
    
    % Initialization of parameters 
    x = ones(k,1);
    y = ones(l,1);
    s = ones(q,1);
    z = ones(q,1);
    
    Z = diag(z);
    S = diag(s);
        
    [m,n] = size(S);
    e = ones(m,1);    
    eta = 0.995;
    
    max_iter = 100  ;
    tol = 10^-9;
    
    % Compute residuals
    rL = H*x + g - A*y - C*z;
    rA = b - A'*x;
    rC = s + d - C'*x;
    rSZ= Z*S*e; 
    mu = z'*s/m;
            
%     disp('-----------')
%     disp(norm(rL,inf))
%     disp(norm(rA,inf))
%     disp(norm(rC,inf))
%     disp(abs(mu))
%     disp('-----------')
    

    Converged = (norm(rL,inf) <= tol) && ...
                (norm(rA,inf) <= tol) && ...
                (norm(rC,inf) <= tol) && ...
                (abs(mu) <= tol);
    i = 1;
%%        
    while ~Converged && (i<max_iter)
%         disp(i)
        % LDL factorization of modified KKT
        s = diag(S);
        z = diag(Z);
        
        % TEMPORARILY FOR GETTING IT TO WORK
        inv_S = inv(S);
        inv_Z = inv(Z);
        
        
        H = H + C*(inv_S*Z)*C';
        KKT = [H, -A; -A', zeros(l,l)];
        [L,D,p] = ldl(KKT,'lower','vector');

        
        %% AFFINE STEP
        % Update residuals
        rL_bar = rL - C * (inv_S*Z) * (rC-inv_Z*rSZ);
        
        % Calculate affine deltas
        tmp = -[rL_bar;rA];
        X(p) = L'\(D\(L\tmp(p)));
        X = reshape(X, length(X), 1);
        dX_aff = X(1:k);
        dY_aff = X(k+1:k+l);
        
        tmp = inv_S*Z;
        dZ_aff = -tmp * C' * dX_aff + tmp * (rC - inv_Z * rSZ);
        dS_aff = -inv_Z * rSZ - inv_Z * S * dZ_aff;
        
        % Compute largest alpha_aff
        idxZ = find(dZ_aff < 0.0);
        idxS = find(dS_aff < 0.0);   
        alpha_aff = min([1.0; -z(idxZ,1)./dZ_aff(idxZ,1); -s(idxS,1)./dS_aff(idxS,1)]);

        % Dual gap
        mu = z'*s/m;
        mu_aff = (z+alpha_aff*dZ_aff)'*(s+alpha_aff*dS_aff)/m;
        ro = (mu_aff/mu)^3;
        
        
        %% CENTERING STEP
        
        % Update residuals
        rSZ_bar= rSZ + dS_aff' * dZ_aff * e - ro * mu * e;
        rL_bar = rL - C * (inv_S*Z) * (rC-inv_Z*rSZ_bar);
        
        % Calculate centering deltas
        tmp = -[rL_bar;rA];
        X = [];
        X(p) = L'\(D\(L\tmp(p)));
        X = reshape(X, length(X), 1);
        dX_cen = X(1:k);
        dY_cen = X(k+1:k+l);
        dZ_cen = - (inv_S*Z) * C' * dX_cen + (inv_S*Z)*(rC - inv_Z * rSZ_bar);
        dS_cen = - inv_Z * rSZ_bar - inv_Z * S * dZ_cen;
        
        
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
       
        Z = diag(z);
        S = diag(s);        
        
        % Update residuals
        rL = H*x + g - A*y - C*z;
        rA = b - A'*x;
        rC = s + d - C'*x;
        rSZ= Z*S*e; 
        mu = z'*s/m;
        
        % Get norms of residuals
%         disp('-----------')
%         disp(norm(rL,inf))
%         disp(norm(rA,inf))
%         disp(norm(rC,inf))
%         disp(abs(mu))
%         disp('-----------')

        Converged = (norm(rL,inf) <= tol) && ...
                    (norm(rA,inf) <= tol) && ...
                    (norm(rC,inf) <= tol) && ...
                    (abs(mu) <= tol);
    
        i = i + 1;
    end

%         disp('-----------')
%         disp(norm(rL,inf))
%         disp(norm(rA,inf))
%         disp(norm(rC,inf))
%         disp(abs(mu))
%         disp('-----------')

        if ~Converged
            x=[];
            y=[];
            z=[];
            s=[];
        end
end