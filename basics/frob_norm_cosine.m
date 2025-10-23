function f = frob_norm_cosine(sig2, p, Kappa_hat)
    

    % Compute first row of Toeplitz matrix from periodic Matérn covariance

R_model =  sig2*toeplitz(cos(2*pi * (0:(p-1)) / p));
    residual = Kappa_hat - R_model;
    f = norm(residual, 'fro'); 

end
