%Compute a feasible starting point x0;
%SetW0 to be a subset of the active constraints at x0;
for k = 0, 1, 2, . . .
%Solve (16.39) to find pk ;
%(16.39) = min_p  ½p^T Gp + g_k^Tp    s.t. a_i^Tp = 0, iW_k
    if pk = 0
    %Compute Lagrange multipliers ˆ?i that satisfy (16.42),
    %with ˆW = Wk ;
    % (16.45) = sum(a_i*lambda_i^ = g = Gx^+ c
        if lamb_i >= 0 for all i %? W_k ? I
            stop %with solution x? = x_k ;
        else
        %j ? arg min_(j?Wk?I) lamb_j ;
        %x_k+1 ? x_k; W_k+1 ? W_k\{j};
    else %(pk ~= 0)
    %Compute ?k from (16.41);
    %(16.41) alfa_k = min(1, min  b_i-a_i^Tx_k/(a_i^Tp_k) 
    %xk+1 ? xk + ?k pk ;
    if %there are blocking constraints
    %ObtainWk+1 by adding one of the blocking
    %constraints toWk ;
    else
    %Wk+1 ? Wk ;
end(for)
