function[x_prior,x_posterior,P_prior,P_posterior,L_prior,L_posterior] = srkf_schmidt(y,x,u,P,A,B,G,C,Q,R,L)
% Square root covariance filter (with QR-decomposition)  from (Schmidt et al., 1968)
    x_posterior = x;
    n = length(x);
    x_prior     = A*x_posterior + B*u;      % prior mean
    M1          = [L'*A'; Q'];
    % Givens rotations
    Upper = upper_tri(M1)
%     [Qq, R]      = qr(M1',0);                % QR decomposition to return upper triangular from the pre-array
%     if size(R,2) > n
%         L_prior     = tril(R',-n);
%     else
%         L_prior = R';
%     end
    L_prior = Upper';
    P_prior     = L_prior*L_prior';                   % prior covariance
    d_y         = y - C*x_prior;            % estimation error  
    M2          = [L_prior'*C'; R'];
    S           = M2'*M2;                   % error covariance
%     [row,col] = find((S<eps*1000));
%     for i=1:(length(row))
%         S(row(i),col(i)) = 0;
%     end
    [L_res,rr]  = chol(S,'lower');
    S_inv       = pinv(S);   
    K           = P_prior*C'*pinv(L_res')*pinv(L_res);         % Kalman gain
    x_posterior = x_prior + K*d_y;          % ostrerior mean  
    L_posterior = (eye(size(L_prior)) - K*C')*L_prior; % L_res*pinv(L_res + R)*
    P_posterior = L_posterior*L_posterior';                   % posterior covariance

    S_new = (abs(2*pi*S));
    den = sqrt(det(S_new));
    num = exp(-0.5*d_y'*S_inv*d_y);
    lik   = num/den;
end