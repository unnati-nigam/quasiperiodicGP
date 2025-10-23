function eqn = dldKappa(n, p,omega,Kappa, X)
       k = floor(n / p); % Number of blocks
    l=n-k*p;
    %p-l
    S = zeros(p, p); % Initialize S matrix
    S_new=zeros(p,p);
    Y = convert_Ydata_vector(X,p);
    Y_matrix=Y(:, 2:end) - omega * Y(:, 1:end-1);
    

    % Loop through blocks to compute S
    for i = 1:(k-1)
        % Extract relevant block segments
        x_current = X((i*p+1):((i+1)*p));
        x_prev = X(((i-1)*p+1):(i*p));
        
        % Compute the residual matrix for this block
        residual = x_current - omega * x_prev;
        
        % Accumulate the outer product
        S = S + residual * residual';
    end
%S
     if mod(n,p)==0
   
        S = S / (k-1);

     else 
         
         S=zeros(p,p);
        %S_new=em(Y_matrix,1000,1e-10);
        Y_matrix1=Y_matrix(1:l,1:(k-1));
        Y_matrix2=Y_matrix((l+1):p,1:(k-1));
            S11=zeros (l,l);
            S12=[];%zeros((p-l),l);
            S21=S12';
            S22=zeros ((p-l), (p-l));

            
       % size(Y_matrix1)
        %size(Y_matrix2)

        S11=(Y_matrix1*Y_matrix1')/(k-1);
        S12=(Y_matrix1*Y_matrix2')/(k-1);
        S21=S12';
        S22=(Y_matrix2*Y_matrix2')/(k-1);
        add=Y_matrix(1:l,k)'*Y_matrix(1:l,k);
        S11star=((k-1)*S11+add)/k;

        S_new(1:l,1:l)=S11star;
        S_new(1:l,(l+1):p)=S11star*inv(S11)*S12;
        S_new((l+1):p,1:l)=(S_new(1:l,(l+1):p))';
        S_new((l+1):p,(l+1):p)=S22-S21*inv(S11)*S12+S21*inv(S11)*S11star*inv(S11)*S12;
       % norm(S_new(1:l,1:l)-S11,'fro')
        %S=S_new
        S=S_new;

    end
    [m, n] = size(S);


% Loop over all diagonals
for k = -(m-1):(n-1)
    % Get linear indices of current diagonal
    [row_idx, col_idx] = find(bsxfun(@minus, (1:m)', 1:n) == k);
    
    % Get values of the current diagonal
    diag_vals = S(sub2ind([m, n], row_idx, col_idx));
    
    % Compute the average of current diagonal
    avg_val = mean(diag_vals);
    
    % Replace all elements on this diagonal with the average
    for i = 1:length(row_idx)
        S(row_idx(i), col_idx(i)) = avg_val;
    end
end

epsilon=1e-6;
f = real(fft(S(1, :)));
f(f<0)=0;
f=abs(f);
f= max(real(f), epsilon) ;%+ 1i * imag(f);
v=ifft(f);
% norm(S/(k-1)-toeplitz(v),'fro')
        
Kappa_est=toeplitz(v);
%eqn=norm(eye(p)-inv(Kappa_est)*S,'Inf');
A=eye(p)-inv(Kappa_est)*S;
eqn=max(A(:));
%norm(S_new-S,'fro')


function Y_matrix = convert_Ydata_vector(Y_data, p)
% Convert long vector Y_data to (p x n) matrix with NaNs for incomplete block
% Inputs:
% - Y_data: long column vector of concatenated Y_i's
% - p: dimension of each full sample
%
% Output:
% - Y_matrix: p x (n+1) matrix with last column padded by NaNs if incomplete

len = length(Y_data);
n_full = floor(len / p);        % number of complete blocks
n_extra = mod(len, p);          % number of entries in the last partial block

% Reshape full part
Y_matrix = reshape(Y_data(1 : n_full*p), p, n_full);

% If there is a leftover partial observation
if n_extra > 0
    last_col = NaN(p, 1);
    last_col(1:n_extra) = Y_data(end - n_extra + 1 : end);
    Y_matrix = [Y_matrix, last_col];  % append as new column
end

end


end
