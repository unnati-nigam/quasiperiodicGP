function pred = pred_element(n, p, w, Kappa, X)
% PRED_ELEMENT Robust element-wise prediction using stable inversion
%   Inputs:
%       n     - total length of series
%       p     - block length
%       w     - AR-like coefficient
%       Kappa - covariance matrix (p x p)
%       X     - input time series (n x 1)
%   Output:
%       pred  - predicted series (1 x n)

pred = zeros(1, n);  % initialize prediction

% Optional: clip scale for numerical stability

for j = 1:n
    if mod(j, p) ~= 0
        i = floor(j / p);
        l = j - i*p;
    else
        i = floor(j / p) - 1;
        l = p;
    end
    
    if j > p
        % AR term
        pred1 = w * X((i-1)*p + l);
        
        % Residual term
        if l > 1
            rhs = X((i*p + 1):(i*p + l-1)) - w * X(((i-1)*p + 1):((i-1)*p + l-1));
            K_sub = Kappa(1:(l-1), 1:(l-1));
            K_row = Kappa(l, 1:(l-1));
            
            % Robust inversion using Cholesky or pseudo-inverse
            try
                L = chol(K_sub, 'lower');       % try Cholesky
                y = L \ rhs;
                x = L' \ y;
            catch
                x = pinv(K_sub) * rhs;         % fallback
            end
            
            pred2 = K_row * x;
        else
            pred2 = 0;
        end
        
    else
        % For j <= p
        if j == 1
            pred1 = X(l);
            pred2 = 0;
        else
            K_sub = Kappa(1:(j-1), 1:(j-1));
            K_row = Kappa(l, 1:(j-1));
            rhs = X(1:(j-1));
            
            try
                L = chol(K_sub, 'lower');
                y = L \ rhs;
                x = L' \ y;
            catch
                x = pinv(K_sub) * rhs;
            end
            
            pred1 = 0;
            pred2 = K_row * x;
        end
    end
    
    % Optional clipping for numerical stability
    pred(j) = pred1 + pred2;
    pred(j) = max(min(pred(j), max(X)), min(X));  % keep within observed data range
end

end

% function pred = pred_element(n, p, w, Kappa, X)
% % Predict the elements based on the input parameters n, p, w, Kappa, and X

% pred = zeros(1, n); % Initialize the prediction vector with zeros, length n

% for j = 1:n  % Loop over each element in the series
%     if mod(j, p) ~= 0  % If j is not a multiple of p
%         i = floor(j / p);  % Determine the block index
%         l = j - i*p;       % Determine the position within the block
%     else  % If j is a multiple of p
%         i = floor(j / p) - 1;  % Adjust block index
%         l = p;                  % Position within block is p
%     end
    
%     if j > p  % If we are beyond the first block
%         pred1 = w * X((i-1)*p + l);  % AR term: w times the previous block element
%         rhs = X((i*p + 1):(i*p + l-1)) - w * X(((i-1)*p + 1):((i-1)*p + l-1));  
%         % Compute the difference between current block segment and AR prediction
        
%         pred2 = Kappa(l, 1:(l-1)) * (Kappa(1:(l-1), 1:(l-1)) \ rhs);  
%         % Solve linear system using Kappa submatrix, then multiply by Kappa row
%     else  % For first block elements (j <= p)
%         if j == 1  % First element of the series
%             pred1 = X(l);  % Prediction is the value itself
%             pred2 = 0;     % No covariance contribution
%         else  % Other elements in first block
%             pred1 = 0;  % No AR contribution
%             pred2 = Kappa(l, 1:(j-1)) * (Kappa(1:(j-1), 1:(j-1)) \ X(1:(j-1)));  
%             % Solve linear system using first j-1 elements
%         end
%     end
    
%     pred(j) = pred1 + pred2;  % Total prediction is sum of AR and covariance contributions
% end
% end

