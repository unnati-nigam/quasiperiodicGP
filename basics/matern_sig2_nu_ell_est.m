function params = matern_sig2_nu_ell_est(p, Kappa_hat)
    % Initial guess: [sig2, nu, ell]
  %initial_guess = [4.8, 2.5, 0.1];%co2
 %lb=[eps,2.5,eps]; ub=[Inf, 2.5, 5];%
 %lb = [eps,eps,eps]; ub=[Inf,1,1]; %co2
    %initial_guess = [4000, 0.01,0.01]; %
    initial_guess = [1, 1.5,2]; %sunspot 
 lb=[eps,1.5,0.1]; ub=[Inf,1.5, 5];% 
 %lb=[eps,eps,eps]; ub=[Inf,Inf, Inf];%sunspot 
%initial_guess = [0.27, 1.5,0.01]; %tide 
 % lb=[eps,1.5,eps]; ub=[Inf, 1.5, 5];
  %lb =[eps,eps,eps];ub=[1,1,1]; %tide 

    % Define the objective function handle
    objective = @(params) frob_norm_matern(params, p, Kappa_hat);

    % Use fminunc or fmincon for constrained/unconstrained optimization
   % options = optimoptions('patternsearch', 'Display','none','AccelerateMesh',true);
options = optimoptions('patternsearch', 'Display','none','AccelerateMesh',true);
params = patternsearch(objective, initial_guess, [], [], [], [], lb, ub, [], options);
end
