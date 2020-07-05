function[x_prior,x_posterior,P_prior,P_posterior] = kf_joseph(y,x,u,P,dfun,mfun,A,C,Q,R)
    x_posterior = x;
    x_prior     = dfun(x_posterior,u);    % prior mean
    P_prior     = A*P*A' + Q;             % prior covariance
    d_y         = y - mfun(x_prior,u);          % estimation error    
    S           = C*P_prior*C' + R;       % error covariance
    S_inv       = pinv(S);
    K           = P_prior*C'*S_inv;       % Kalman gain
    x_posterior = x_prior + K*d_y;        % ostrerior mean          
    P_posterior = (eye(length(x)) - K*C)*P_prior*(eye(length(x)) - K*C)' + K*R*K'; % posterior covariance
        
%     S_new = (abs(2*pi*S));
%     den = sqrt(det(S_new));
%     num = exp(-0.5*d_y'*S_inv*d_y);
%     lik   = num/den;
end