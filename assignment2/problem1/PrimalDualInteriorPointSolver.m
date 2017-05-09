function [x,mu,lambda,Converged] = PrimalDualInteriorPointSolver(g,A,b,x)
    [k,l] = size(x);
    [m,n] = size(A);

    
    % Initialization of parameters
    lambda  = ones(n,1);
    mu      = zeros(m,1);
    eta     = 0.995;
    max_iter= 100;
    tol     = 10^-8;
    
    
    % Compute residuals
    rA = A*x-b;   
    rC = x'*lambda;
    rL = g-A'*mu-lambda;                             
    s  = rC/k;              


    Converged = (norm(rA,inf) <= tol) && ...
                (norm(rL,inf) <= tol) && ...
                (abs(s) <= tol);
    
    %%
    i = 1;
    while ~Converged && (i<max_iter)
        % Cholesky factorization of normal equations
        G = A * diag(x./lambda) * A';        
        [L,p,r] = chol(sparse(G),'lower','vector');
                
        
        %% AFFINE STEP
        
        % Calculate affine deltas
        tmp_1 = (x.*rL + rC)./lambda;
        tmp_2 = A*tmp_1 - rA;
        dMu(r) = L'\(L\tmp_2(r));
        dMu = reshape(dMu, length(dMu),1);
        dX_aff = x./lambda.*(A'*dMu) - tmp_1;
        dLambda = -(rC+lambda.*dX_aff)./x;

        
        % Compute largest alpha_aff
        idxX      = find(dX_aff < 0.0);
        idxLambda = find(dLambda < 0.0);
        alpha_aff = min([1.0; -x(idxX,1)./dX_aff(idxX,1)]);
        beta_aff  = min([1.0; -lambda(idxLambda,1)./dLambda(idxLambda,1)]);
        
        
        x_aff = x + alpha_aff*dX_aff;
        lambda_aff = lambda + beta_aff*dLambda;

        
        % Dual gap
        s_aff = x_aff' * lambda_aff / n;
        sigma = (s_aff / s)^3;    
        tau   = sigma * s;    

        
        %% CENTERING STEP
        
        % Update residual
        rC = rC + dX_aff.*dLambda - tau;
        
        
        % Calculate centering deltas
        tmp_1  = (x.*rL + rC)./lambda;
        tmp_2  = -rA + A*tmp_1;
        dMu    = L'\(L\tmp_2);
        dX_cen = x./lambda.*(A'*dMu) - tmp_1;
        dLambda= -(rC+lambda.*dX_cen)./x;

        
        % Compute largest alpha_cen
        idxX      = find(dX_cen < 0.0);
        idxLambda = find(dLambda < 0.0);
        alpha_cen = min([1.0; -x(idxX,1)./dX_cen(idxX,1)]);
        beta_cen  = min([1.0; -lambda(idxLambda,1)./dLambda(idxLambda,1)]);

        
        % Update position
        x = x + (eta*alpha_cen)*dX_cen;
        mu = mu + (eta*beta_cen)*dMu;
        lambda = lambda + (eta*beta_cen)*dLambda;

        
        % Update residuals   
        rA = A*x - b;               
        rC = x.*lambda;   
        rL = g - A'*mu - lambda; 
        s  = x'*lambda/k;

        Converged = (norm(rA,inf) <= tol) && ...
                    (norm(rL,inf) <= tol) && ...
                    (abs(s) <= tol);
    
        i = i+1;
        
    end

    if ~Converged
        x = [];
        mu = [];
        lambda = [];
    end
end
    