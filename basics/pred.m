function pred = pred(n, p, w, X)
    % Predict the elements based on the input parameters n, p, w, R, and X
    pred = zeros(1, n); % Initialize the prediction vector
    for j = 1:n 
         if mod(j, p) ~= 0
            i = floor(j / p);
            l = j-i*p;
        else
            i = floor(j / p) - 1;
            l = p;
        end
        if j>p
            pred(j) = w * X((i-1)*p + l);
        else 
            pred(j) = 0;
        end
    end 
end 
