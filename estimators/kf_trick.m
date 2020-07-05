function[x_prior,x_posterior,P_prior,P_posterior] = kf_trick(y,x,u,P,A,B,C,Q,R)
    x_posterior = x;
    x_prior     = A*x_posterior + B*u;    % prior mean
    P_prior     = A*P*A' + Q;             % prior covariance
    d_y         = y - C*x_prior;          % estimation error    
    S           = C*P_prior*C' + R;       % error covariance
    S_inv       = pinv(S);
    K           = P_prior*C'*S_inv;       % Kalman gain
    x_posterior = x_prior + K*d_y;        % ostrerior mean          
    P_posterior = (eye(4) - K*C)*P_prior; % posterior covariance
    P_posterior = (P_posterior + P_posterior')/2;
      
    S_new = (abs(2*pi*S));
    den = sqrt(det(S_new));
    num = exp(-0.5*d_y'*S_inv*d_y);
    lik   = num/den;
end