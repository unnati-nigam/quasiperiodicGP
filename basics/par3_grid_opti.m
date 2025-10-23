function result = par3_grid_opti(n, p, X)

% Define grid
omega_grid = 0:0.01:0.99;
theta_grid = 0.5:0.01:1.5;
sig2_grid  = 0.5:0.01:1.5;

% Preallocate
min_logL = Inf;
best_params = [];

% Sin term (only needs to be computed once)
s = sin(pi * (0:p-1) / p);

% === Grid Search ===
for i = 1:length(omega_grid)
    omega = omega_grid(i);

    for j = 1:length(sig2_grid)
        sig2 = sig2_grid(j);

        for k = 1:length(theta_grid)
            theta = theta_grid(k);

            % Build Kappa
            kappa_vec = sig2 * exp(-theta^2 * s.^2);
            Kappa = toeplitz(kappa_vec);

            % Evaluate log-likelihood
            try
                logL = logL_exact(n, p, omega, Kappa, X);
            catch
                logL = Inf;  % Safeguard against numerical failure
            end

            % Track best result
            if isfinite(logL) && logL < min_logL
                min_logL = logL;
                best_params = [omega, sig2, theta];
            end
        end
    end
end


omega_best = best_params(1);
sig2_best  = best_params(2);
theta_best = best_params(3);

result=[omega_best, theta_best,sig2_best];

