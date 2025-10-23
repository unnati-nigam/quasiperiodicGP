function obj = log_L_TSP(n,para, Y,corr,regr,P) % using smt package
    delta = para(1);     theta = para(2);     omega = para(3); 
    data = struct('corr',corr, 'regr',regr, 'P', P, 'Y',Y);
    
        p = P;    k = floor(n/p);    p1 = n - p*k;
        Y = zeros(p,k);    Y(1:k*p) = data.Y(1:k*p);    Ys = data.Y(k*p+1:n);
        [So, Uo] = kmseig(k, omega);
        r0 = data.corr(delta, theta, p, (1:p)', 1);
        r0(1) = r0(1) + eps;    Sr = real(fft(r0));    Sr(Sr<0) = 0; 
        S = Sr * So' +delta^2;   
        Lambda = 1./S;
        
        UrYUo = fft(Y*Uo)/sqrt(p);
        ChiYY = real(sum(sum(conj(UrYUo).*(UrYUo.*Lambda))));
       
        if p1 == 0
            sigma2 =ChiYY/n;
            likelihood = (n*log(sigma2) + sum(sum(log(real(S)))) + n + n*log(2*pi))/2;
        else
             r = r0;             r(1) = r(1) + delta^2; 
            UoW = Uo'*omega.^(k:-1:1)';
            r = r - ifft((Lambda*(UoW.^2)).*(Sr.^2));
            Yd = ifft(Sr.*((UrYUo.*Lambda)*UoW))*sqrt(p);
            Yd = real(Ys - Yd(1:p1)); 
            L = toep_chol(r(1:p1)); 
            Ydl = L \ Yd;   
            YPiY = Ydl'*Ydl;
            sigma2 = (ChiYY+YPiY)/n;
            likelihood = (n*log(sigma2) + sum(sum(log(real(S)))) + sum(2*log(diag(L))) + n + n*log(2*pi))/2;
                
        end
          obj = likelihood;
end