function f = frob_norm_matern(sig2_nu_ell, p, Kappa_hat)
    % Unpack parameters
    sig2 = sig2_nu_ell(1);
    nu   = sig2_nu_ell(2);
    ell  = sig2_nu_ell(3);

    % Compute first row of Toeplitz matrix from periodic Matérn covariance
    d = 0:(p-1);
    row = sig2 * matern_cov(d, nu, ell, p);  % matern_cov must be defined

    % Construct Toeplitz matrix
    R = toeplitz(row);

    % Set diagonal to sig2 explicitly (as in your R code)
    R(1:p+1:end) = sig2;

    % Compute Frobenius norm of the difference
    f = sqrt(sum((Kappa_hat(:) - R(:)).^2));
end
