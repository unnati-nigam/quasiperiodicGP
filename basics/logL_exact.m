function logL = logL_exact(n, p, omega, Kappa, X)
 
   % [U, E] = eig(Kappa);                % U: eigenvectors, E: eigenvalues (diagonal)
   % eig_values = real(diag(E));     % Extract real part of eigenvalues
   % Kappa_inv = real(U * diag(1 ./ eig_values) * U'); % Compute R_inv using eigen-decomposition
    %R_inv = real(U * diag(1 ./ eig_values) * U'); % Compute R_inv using eigen-decomposition
    
    epsilon = 1e-10;
    Kappa_reg = Kappa + epsilon * eye(p);
    L = chol(Kappa_reg, 'lower');
    Kappa_inv = L' \ (L \ eye(p));               % Equivalent to inv(R)
    log_detKappa = 2 * sum(log(diag(L)));
    % Initialize sum
    sum_term = 0;

    % Compute the number of blocks
    k = floor(n / p);

    % Loop over the blocks
    for i = 2:k
        % Extract relevant parts of X
        X_curr = X(((i-1)*p+1):(i*p));
        X_prev = X(((i-2)*p+1):((i-1)*p));
        diff = X_curr - omega * X_prev;
       % i

        % Add to the sum
        sum_term = sum_term + (diff') * Kappa_inv * diff;
    end
    sum_term1=(X(1:p))'* Kappa_inv*X(1:p);
  %  detKappa = prod(eig_values);
        
   logL = (n*log(2*pi))/2 +(k*log_detKappa)/2-p*log(1-omega^2)/2 + sum_term/2+(1-omega^2)*sum_term1/2;
    end
