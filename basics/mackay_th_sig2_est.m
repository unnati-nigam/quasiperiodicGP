function params = mackay_th_sig2_est(p, R)
    % Estimates parameters theta and sigma^2 using optimization
    % Inputs:
    %   n - Number of observations
    %   p - Block size
    %   R - Covariance matrix
    %   X - Data vector
    % Output:
    %   params - Optimized parameters [theta, sigma^2]

 %initial_guess = [2.15,4.8];  %co2
% lb=[eps, eps]; ub=[Inf, Inf];%  
 %lb=[1, 3]; ub=[4,10]; %co2 
    %initial_guess = [12, 0.27];%tide
  %lb=[eps,eps]; ub=[Inf, Inf];    %
 %lb=[eps,0.25]; ub=[Inf, 1]; %tide1
  % initial_guess=[4598,0.01];% 
  initial_guess=[12,0.27];%sunspot(month)[4000,0.01];
    lb=[eps,eps];ub=[Inf,Inf]; %sunspot
   %initial_guess=[1,1]; %simulation
   %lb=[eps,eps]; ub=[Inf, Inf];%sim
    objective = @(params) frob_norm_mackay(params, p, R);
    
    
    % Perform optimization
%    params = fmincon(objective, initial_guess, [], [], [], [], lb,ub,[], options);
options = optimoptions('patternsearch', 'Display','none','AccelerateMesh',true);
params = patternsearch(objective, initial_guess, [], [], [], [], lb, ub, [], options);

end
