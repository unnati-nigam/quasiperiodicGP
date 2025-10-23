function result = par_est(n, p_search, X)
    % par_est estimates parameters based on log-likelihood maximization
    % Inputs:
    %   n        - Data length
    %   p_search - Vector of candidate values for parameter p
    %   X        - Data matrix
    % Outputs:
    %   result   - A structure containing:
    %              p_est: Estimated p
    %              w_cap: Estimated w
    %              R_cap: Estimated R

    % Initialize storage
    l = zeros(1, length(p_search)); % Log-likelihood values

    % Loop over the candidate values for p
    for i = 1:length(p_search)
        % Compute w and R using w_R_est function
        var = w_Kappa_est(n, p_search(i), X);
        w = var.w;
        Kappa = var.Kappa;
       
        % Compute log-likelihood
       % p_search(i)
        l(i) = logL(n, p_search(i), w, Kappa, X);
        %logL(n, p_search(i), w, Kappa, X)
    end

    % Find the index of the minimum log-likelihood
    [~, ind] = min(l);

    % Get the best p value
    p_est = p_search(ind);

    % Compute w and R for the estimated p
    var = w_Kappa_est(n, p_est, X);
    w_cap = var.w;
    Kappa_cap = var.Kappa;
%    sig2_cap=var.sig2;

    % Return the results in a structure
    result = struct('p', p_est, 'w', w_cap, 'Kappa', Kappa_cap);
end
