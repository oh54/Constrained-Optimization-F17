function [x, lambda] = NullSpaceQR(n, u_bar, d_0)
    [H,g,A,b] = ConstructFactorization(n, u_bar, d_0);
    n = size(A,1);
    m = size(A,2); 
    
    [Q,R_bar] = qr(A);
    Q_1 = Q(:,1:m); 
    Q_2 = Q(:,m+1:n); 
    R = R_bar(1:m,1:m);
    
    X_y = R'\b;
    X_z = (Q_2'*H*Q_2) \ (-Q_2'*(H*Q_1*X_y+g));
    x = Q_1*X_y + Q_2*X_z;
    lambda = R .\ (Q_2'*(H*x+g));
end
