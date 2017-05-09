function [ mu ] = get_mu( mu, grad_fk, p_k, B_k, c_k, sigma, rho, margin )
    num = grad_fk' * p_k + (sigma/2) * p_k' * B_k * p_k;
    den = (1 - rho) * norm(min(0, c_k), 1);
    mu_candidate = num/den;

    if mu < mu_candidate
        mu = mu_candidate + margin;
    end
end

