function eqn = dldw(n, p,w,R, X)
    k = floor(n / p); % Number of blocks

    % Spectral decomposition of R
    [U, D] = eig(R); 
    e = diag(D); % Eigenvalues
    R_inv = U * diag(1 ./ e) * U'; % Inverse of R using spectral decomposition

    % Compute numerator (i1i) and denominator (ii)
    i1i = 0; % Initialize numerator
    for i = 2:k
        x1 = X(((i-2)*p+1):((i-1)*p));
        x2 = X(((i-1)*p+1):(i*p));
        i1i = i1i + x1' * R_inv * x2;
    end

    ii = 0; % Initialize denominator
    for i = 1:(k-1)
        x = X(((i-1)*p+1):(i*p));
        ii = ii + x' * R_inv * x;
    end
%i1i/ii
        Rstar=R(1:(n-k*p),1:(n-k*p));
        [U, D] = eig(Rstar); 
        e = diag(D); % Eigenvalues
        Rstar_inv = U * diag(1 ./ e) * U'; % Inverse of R using spectral decomposition
        xi=X((k*p-p+1):(n-p));
        xi1=X((k*p+1):n);
        nstar=xi'*Rstar_inv*xi1;
        dstar=xi'*Rstar_inv*xi;
    % Compute w estimate
    if mod(n,p)==0
        nstar=0; dstar=0;
    end
eqn=(ii+dstar)*w-(i1i+nstar);
end
