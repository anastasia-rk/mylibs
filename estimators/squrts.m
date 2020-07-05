function[x_s,p_s,lik_s] = squrts(x_sm,x,u,P_sm,P_p,Q,dfun)
%% preliminaries
L       = length(x);  % state vector dimension
alpha   = 0.0001;       % scaling constant - from literature (set to larger value for wider distribution spread)
kappa   = 0;      % from literature
beta    = 2;          % optimal in case distribution of x is Gaussian
lambda  = alpha^2*(L + kappa) - L;  % scaling factor for sigma points
c       = L + lambda;               % scaling factor for covariance
x_prior = zeros(L,1);
P_prior = zeros(L,L);
C_xx    = zeros(L,L);

% State augmentation
Qq = Q; %(3:4,3:4);
m1 = size(P_p,1);
m2 = size(Qq,1);
O = zeros(m1,m2);
P = [P_p O; O' Qq];
qq = zeros(size(Qq,1),1);
x_tilde = [x' qq']';
%% Sigma points
[LM,DM,E,pneg] = mchol(P);
S_tilde = sqrt(c)*LM*sqrt(DM);
% first values
X_tilde{1} = x_tilde;
K = numel(x_tilde);
X{1} = x;
q{1} = qq;
% placeholders for sigmas
xxx_t = x_tilde(:,ones(1,numel(x_tilde)));
X3 = [x_tilde xxx_t+S_tilde xxx_t-S_tilde]; 
for i = 2:2*K+1
    X_tilde{i} = X3(:,i);
    X{i} = X_tilde{i}(1:L,:);
    q{i} = X_tilde{i}(L+1:end,:);
end
% Weights
W_s(1) = lambda/c;
W_c(1) = lambda/c + (1 - alpha^2 + beta);
W_s (2:2*K+1) = 1/(2*c);  % for states
W_c (2:2*K+1) = 1/(2*c);  % for covariances


%% Unscented transform
for i=1:2*K+1
    Xx{i}     = dfun(X{i},u) + q{i}; % sigma state propagation
    x_prior   = x_prior + W_s(i)*Xx{i};                                    % prior estimate
end

for i=1:2*K+1
    P_prior   = P_prior + W_c(i)*(Xx{i} - x_prior)*(Xx{i} - x_prior)'; % prior covariance
    C_xx      = C_xx + W_c(i)*(X{i} - x)*(Xx{i} - x_prior)';           % cross-covariance
end

%% Correction
S_inv = pinv(P_prior);
d_x = x_sm - x_prior;

D    = C_xx*pinv(P_prior);   % Kalman gain
x_s = x + D*d_x;          % State correction
p_s = P_p - D*(P_sm - P_prior)*D';   % Covariance correction


% Calculate likelihood for IMM

S_new = (abs(2*pi*P_prior));
den = sqrt(det(S_new));
num = exp(-0.5*d_x'*pinv(P_prior)*d_x);
lik_s   = num/den;