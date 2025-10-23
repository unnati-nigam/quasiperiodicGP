function pred = actual_prediction(n, p, w, Kappa, y)
    k=floor(n/p);
    [i, j] = ndgrid(1:k, 1:k);
    Omega = w.^(abs(i - j));
    
    Sigma=kron(Omega, Kappa)/(1-w^2);
    pred=zeros(n,1);
    pred(1) = y(1);
    for t=2:n
        pred(t)=Sigma(1:t-1, t)' * (Sigma(1:t-1, 1:t-1) \ y(1:t-1));
    end
end
