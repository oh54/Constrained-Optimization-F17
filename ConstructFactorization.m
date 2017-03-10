function [H,g,A,b] = ConstructFactorization(n, u_bar, d_0)

    % Create diagonal 1 and -1 above the main diagonal
    tmp_1 = [diag(ones(n-1,1)), zeros(n-1,1)];
    tmp_2 = [zeros(n-1,1), diag(-ones(n-1,1))];
    A = tmp_1 + tmp_2;

    % Get rid of paddings necessary for matrices consistency
    A(end,:) = [];
    A(:,end) = [];

    % Append 2 zero's columns to the end
    A(:,n:n+1) = zeros(n-2,2);

    % Add the first constraint
    tmp_3 = zeros(1,n+1);
    tmp_3([1,n]) = [-1, 1];
    A = [tmp_3; A];
    
    % Add last constraint  
    tmp_4 = zeros(1,n-1);
    tmp_4(n-1:n+1) = [1, -1, -1];
    A = [A; tmp_4]';
    b = zeros(n, 1);
    b(1) = -d_0;
    H = diag(ones(1,n+1));
    g = ones(n+1,1)*2*u_bar;

end