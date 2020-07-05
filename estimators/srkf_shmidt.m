function[x_prior,x_posterior,P_prior,P_posterior,L_prior,L_posterior] = srkf_schmidt(y,x,u,P,A,B,C,Q,R,L)
% Square root covariance filter (with QR-decomposition)  from (Schmidt et al., 1968)
    x_posterior = x;
    x_prior     = A*x_posterior + B*u;      % prior mean
    M1          = [L'*A'; Q'];
    [Q, R]      = qr(M1',0);                % QR decomposition to return upper triangular
    L_prior     = R';
    P_prior     = L_prior*L_prior';                   % prior covariance
    d_y         = y - C*x_prior;            % estimation error  
    M2          = [L_prior'*C'; R'];
    S           = M2*M2';                   % error covariance
    if all(eig(S)> eps)
        L_res       = chol(S,'lower'); 
        S_inv       = pinv(S);
        K           = P_prior*C'*S_inv;         % Kalman gain
        x_posterior = x_prior + K*d_y;          % ostrerior mean  
        L_posterior = (eye(size(S1)) - K_t*L_res*pinv(L_res + R)*C')*L_prior;
        P_posterior = L_posterior*L_posterior';                   % posterior covariance
    else
    end
    S_new = (abs(2*pi*S));
    den = sqrt(det(S_new));
    num = exp(-0.5*d_y'*S_inv*d_y);
    lik   = num/den;
    
end