function [w_boot, Kappa_general_boot, Kappa_mackay_boot,Kappa_matern_boot, Kappa_cosine_boot] = bootstrap_par_se(y, p)
%BOOTSTRAP_SE Compute bootstrap standard errors and 95% CI of w and first row of Kappa
%
%   INPUTS:
%       y - time series vector (n x 1)
%       p - block length / period
%
%   OUTPUTS:
%       SE_w         - standard error of w (1 x length(w))
%       SE_Kappa     - standard errors of the first row of Kappa (1 x p)
%       w_lower      - 2.5% quantile of w (1 x length(w))
%       w_upper      - 97.5% quantile of w (1 x length(w))
%       Kappa_lower  - 2.5% quantile of first row of Kappa (1 x p)
%       Kappa_upper  - 97.5% quantile of first row of Kappa (1 x p)

B = 1000;      
rng(1);
n = length(y);
k = floor(n/p);

% 1. Estimate parameters from original data
par = par_est(n, p, y);
w_hat = par.w;

% 2. Compute residuals block-wise
residuals = zeros(k-1, p);
for i = 1:(k-1)
    y_prev = y((i-1)*p + 1 : i*p);
    y_curr = y(i*p + 1 : (i+1)*p);
    residuals(i,:) = y_curr(:)' - w_hat * y_prev(:)';
end

residuals = residuals - mean(residuals, 'all'); % center residuals
%residuals = residuals * sqrt(var(y(:)) / var(residuals(:)));  % rescale

% Initialize storage
w_boot = zeros(B, length(w_hat));
Kappa_general_boot = zeros(p,p, B);
Kappa_mackay_boot = zeros(2,B);
Kappa_matern_boot = zeros(3,B);
Kappa_cosine_boot = zeros(1,B);

% 3. Bootstrap resampling
for b = 1:B
    idx = randi(k-1, k-1, 1);  
    z_star = residuals(idx, :);
    
    % Construct bootstrap series
    yb = zeros(1, n);
    yb(1:p) = y(1:p);
    for i = 1:(k-1)
        yb((i*p + 1):((i+1)*p)) = w_hat * yb(((i-1)*p + 1):(i*p)) + z_star(i,:);
    end
    b
    % Re-estimate parameters from bootstrap sample
    par_b = par_est(n, p, yb');
    w_boot(b,:)       = par_b.w;
    Kappa_general_boot(:,:,b) = par_b.Kappa;
    Kappa_mackay_boot(:,b) = mackay_th_sig2_est(par.p, par_b.Kappa);
    Kappa_matern_boot(:,b) = matern_sig2_nu_ell_est(par.p, par_b.Kappa);
    Kappa_cosine_boot(:,b) = cosine_sig2_est(par.p, par_b.Kappa);
end

% 4. Compute bootstrap standard errors
%SE_w = std(W_boot, 0, 1);    
%SE_Kappa = squeeze(std(Kappa_boot(1,:,:), 0, 3))';  % 1 x p

% 5. Compute 2.5% and 97.5% quantiles using quantile
%w_lower = quantile(W_boot, 0.025, 1);
%w_upper = quantile(W_boot, 0.975, 1);

%Kappa_lower = squeeze(quantile(Kappa_boot(1,:,:), 0.025, 3))';  % 1 x p
%Kappa_upper = squeeze(quantile(Kappa_boot(1,:,:), 0.975, 3))';  % 1 x p

end


