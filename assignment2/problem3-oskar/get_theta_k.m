function [ val ] = get_theta_k( B_k, s_k, y_k )
    
    % quadratic term of the approximation
    quad = s_k' * B_k * s_k;

    if s_k' * y_k >= 0.2 * quad
        val = 1;
    else
        val = (0.8 * quad) / (quad - s_k' * y_k);
    end
end

