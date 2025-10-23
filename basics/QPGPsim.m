function X = QPGPsim(n, p, theta, omega, sig2)
    % QPGPsim generates a quasi-periodic Gaussian process simulation
    % Inputs:
    %   n     - Length of the desired output
    %   p     - Periodicity parameter
    %   theta - Parameter controlling the covariance structure
    %   omega - Weighting factor for temporal dependence
    %   sig2  - Variance parameter
    % Output:
    %   X     - Simulated process

    % Construct the Toeplitz covariance matrix
    Kappa = sig2 * toeplitz(exp(-theta^2 * (sin(pi * (0:(p-1)) / p)).^2));

    % Number of segments to generate
    k = floor(n / p);

    % Preallocate the output vector
    if mod(n,p)==0
    
    X = zeros(n, 1);

    % Initial segment
    X(1:p) = mvnrnd(zeros(1, p), Kappa / (1 - omega^2));
    % Iteratively generate segments
    for i = 2:k
        X(((i-1)*p+1):(i*p)) = omega * X(((i-2)*p+1):((i-1)*p))' + mvnrnd(zeros(1, p), Kappa);
    end

    else 
    X=zeros((k+1)*p,1);
    X(1:p) = mvnrnd(zeros(1, p), Kappa / (1 - omega^2));
    for i = 2:(k+1)
        X(((i-1)*p+1):(i*p)) = omega * X(((i-2)*p+1):((i-1)*p))' + mvnrnd(zeros(1, p), Kappa);
    end
    X=X(1:n);
   
end
