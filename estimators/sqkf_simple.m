function[x_prior,x_posterior,P_prior,P_posterior,L_prior,L_posterior] = srkf_simple(y,x,u,P,dfun,mfun,A,C,Q,R,L)
% Square root KF simple form (without ortogonal rotation) from (Carraro&Sartore, 1987)
    Q           = chol(Q,'lower');
    R           = chol(R,'lower');
    x_posterior = x;
    x_prior     = dfun(x_posterior,u);      % prior mean
    M1          = [L'*A'; Q'];
    P_prior     = M1'*M1;                   % prior covariance
    [L1,D1,E,pneg] = mchol(P_prior);
%     L_prior     = chol(P_prior,'lower');    % note that matlab outputs upper triangule
    L_prior     = L1*sqrt(D1);
    d_y         = y - mfun(x_prior,u);            % estimation error  
    M2          = [L_prior'*C'; R'];
    S           = M2'*M2;                   % error covariance
    S_inv       = pinv(S);
    K           = P_prior*C'*S_inv;         % Kalman gain
    x_posterior = x_prior + K*d_y;          % ostrerior mean  
    M3          = [L_prior'*(eye(4) - K*C); R'*K'];
    P_posterior = M3'*M3;                   % posterior covariance
    [L2,D2,E,pneg] = mchol(P_posterior);
    L_posterior = L2*sqrt(D2);
%     L_posterior = chol(P_posterior,'lower');                  % note that matlab outputs upper triangule
%     S_new = (abs(2*pi*S));
%     den = sqrt(det(S_new));
%     num = exp(-0.5*d_y'*S_inv*d_y);
%     lik   = num/den;
end