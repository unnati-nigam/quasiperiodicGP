function matern_cov = matern_cov(d, nu, ell, p)
    % MATLAB equivalent of the R function matern_cov

    if nu == 0.5
        matern_cov = exp(-sqrt((sin(2*pi*d/p)).^2) / ell);
    else
        factor1 = (2^(1 - nu)) / gamma(nu);
        factor2 = ((sqrt(2 * nu) * 2*sqrt((sin(pi*d/p)).^2)) / ell).^nu;
        besselk_val = besselk(nu, (sqrt(2 * nu) * 2*sqrt((sin(pi*d/p)).^2)) / ell);
        matern_cov = factor1 * factor2 .* besselk_val;
    end
end
