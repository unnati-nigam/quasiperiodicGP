function params = cosine_sig2_est(p, R)
    % Estimates parameters theta and sigma^2 using optimization
    % Inputs:
    %   n - Number of observations
    %   p - Block size
    %   R - Covariance matrix
    %   X - Data vector
    % Output:
    %   params - Optimized parameters [theta, sigma^2]
initial_guess=4000;%0.5;% sunspot
lb=eps; ub=Inf;%sunspot
%initial_guess=5;% co2
%lb=eps; ub=Inf;%
%initial_guess = 0.27; %tide;
%lb=eps; ub=Inf;%tide
    objective = @(params) frob_norm_cosine(params, p, R);
    
    
    % Perform optimization
%    params = fmincon(objective, initial_guess, [], [], [], [], lb,ub,[], options);
options = optimoptions('patternsearch', 'Display','none','AccelerateMesh',true);
params = patternsearch(objective, initial_guess, [], [], [], [], lb, ub, [], options);


end
