function [ val ] = q_mu( f_k, grad_fk, hess_Lk, mu, p, c_k, A_k )


    m_p = norm(min(0, c_k + A_k * p), 1);
    val = f_k + grad_fk' * p + 1/2 * p' * hess_Lk * p + mu * m_p;
    

end

