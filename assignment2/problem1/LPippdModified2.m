function [x,info,lambda,s,iter] = LPippdModified(c,A,b,x)

[m,n]=size(A);

maxit = 100;
tolL = 1.0e-9;
tolA = 1.0e-9;
tols = 1.0e-9;

eta = 0.99;

s = ones(n,1);
lambda = zeros(m,1);

% Compute residuals
r_c = A'*lambda + s - c;    % Lagrangian gradient
r_b = A*x - b;               % Equality Constraint
xs = x.*s;             % Complementarity
mu = sum(xs)/n;              % Duality gap (measure?)

% Converged
Converged = (norm(r_c,inf) <= tolL) && ...
            (norm(r_b,inf) <= tolA) && ...
            (abs(mu) <= tols);

%
%x_k = x;
%lambda_k = lambda;
%s_k = s;

iter = 0;
while ~Converged && (iter<maxit)
    iter = iter+1;
    
    %x = x_k;
    %lambda = lambda_k;
    %s = s_k;
    
    %X = diag(x);
    %I = eye(size(X));
    %S = diag(s);
    %e = ones(n);
    
    % 14.30
    %F = [ -r_c ; -r_b ; -X*S*e ];
    %J_F = [ 0 A' I ; A 0 0 ; S 0 X ];
    % 2*n + m sized vector
    %affs = J_F \ F;
    
    % calculate alpha_aff_pri and alpha_aff_dual
    % calculate mu_aff
    % set sigma
    % solve for delta_x delta_lambda delta_s
    % calculate alpha_pri and alpha_dual
    % update x,lambda,s
    
    % ====================================================================
    % Form and Factorize Hessian Matrix
    % ====================================================================
    xdivs = x./s;
    H = A*diag(xdivs)*A';
    L = chol(H,'lower');
    
    % ====================================================================
    % Affine Step
    % ====================================================================
    % Solve
    tmp = (x.*r_c + xs)./s;
    rhs = -r_b + A*tmp;
    
    ds = L'\(L\rhs);
    dx = xdivs.*(A'*ds) - tmp;
    ds = -(xs+s.*dx)./x;
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(ds < 0.0);
    beta = min([1.0; -s(idx,1)./ds(idx,1)]);
    
    % ====================================================================
    % Center Parameter
    % ====================================================================
    xAff = x + alpha*dx;
    sAff = s + beta*ds;
    sAff = sum(xAff.*sAff)/n;
    
    sigma = (sAff/mu)^3;
    tau = sigma*mu;    

    % ====================================================================
    % Center-Corrector Step
    % ====================================================================
    xs = xs + dx.*ds - tau;
    
    tmp = (x.*r_c + xs)./s;
    rhs = -r_b + A*tmp;
    
    ds = L'\(L\rhs);
    dx = xdivs.*(A'*ds) - tmp;
    ds = -(xs+s.*dx)./x;
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(ds < 0.0);
    beta = min([1.0; -s(idx,1)./ds(idx,1)]);

    % ====================================================================
    % Take step 
    % ====================================================================
    x = x + (eta*alpha)*dx;
    lambda = lambda + (eta*beta)*ds;
    s = s + (eta*beta)*ds;
    
    % ====================================================================
    % Residuals and Convergence
    % ====================================================================
    % Compute residuals
    r_c = c - A'*lambda - s;    % Lagrangian gradient
    r_b = A*x - b;               % Equality Constraint
    xs = x.*s;             % Complementarity
    mu = sum(xs)/n;              % Duality gap

    % Converged
    Converged = (norm(r_c,inf) <= tolL) && ...
                (norm(r_b,inf) <= tolA) && ...
                (abs(mu) <= tols);
end

%%
% Return solution
info = Converged;
if ~Converged
    x=[];
    lambda=[];
    s=[];
end
    