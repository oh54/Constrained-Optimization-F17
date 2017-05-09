function [ decrease_term ] = get_decrease_term( mu, eta, alpha_k, grad_fk, p_k, c_k  )
% p.540 18.28 & p.541 18.29 
    dir_grad = grad_fk' * p_k - mu * norm(min(0, c_k), 1);
    decrease_term = eta * alpha_k * dir_grad;
end

