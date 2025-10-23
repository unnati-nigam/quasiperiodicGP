function var_out = pred_var(n, p, w, Kappa)
    var_out = zeros(n,1);  % Preallocate
    
    % Precompute inverses of all leading principal submatrices
    K_inv = cell(p,1);
    for l = 2:p
        K_inv{l} = Kappa(1:(l-1), 1:(l-1)) \ eye(l-1);  % Avoid inv()
    end
    
    for j = 1:n
        if mod(j, p) ~= 0
            i = floor(j / p);
            l = j - i*p;
        else
            i = floor(j / p) - 1;
            l = p;
        end
        
        if l > 1
            Kappa_lt1 = Kappa(l, 1:(l-1));  % Row
            Kappa_lt_inv = K_inv{l};
            cross_term = Kappa_lt1 * Kappa_lt_inv * Kappa_lt1';
        else
            cross_term = 0;
        end
        
        if j > p
            v = (w^2) * Kappa(l,l) / (1 - w^2) + (1 - w^2 + 2*w) * cross_term;
        else
            v = cross_term / (1 - w^2);
        end
        
        var_out(j) = v;
    end
end
