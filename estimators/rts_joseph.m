function[x_s,P_s] = rts_joseph(x_sm,x_pr,x_p,P_sm,P_pr,P_p,u_p,A,B,C,Q,R)
    x_t1T      = x_sm; 
    x_tt       = x_p;  % posterior estimate from Kalman filter
    x_t1t      = A*x_p + B*u_p; % prior estimate from Kalman filter
    P_t1T(:,:) = P_sm;
    P_tt (:,:) = P_p;  % posterior covariance from KF
    P_t1t(:,:) = A*P_p*A' + Q; %  prior covariance from KF
    S_t = P_tt*A'*pinv(P_t1t); % tb denotes time t=t-1
    x_s = x_tt + S_t*(x_t1T - x_t1t); % + P_tt*C'*pinv(R + C*P_tt*C)
    P_s = P_tt  +  S_t*(P_t1T - P_t1t)*S_t'; % P_tt*C'*pinv(R + C*P_tt*C')*C*P_tt'
     

