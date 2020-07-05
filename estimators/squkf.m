function[x_posterior,S_posterior] = squkf(y,x,u,S,Q,R,dfun,mfun)
% preliminaries
L       = length(x);  % state vector dimension
M       = length(y);  % measurement vector dimension
alpha   = 0.01;       % scaling constant - from literature (set to larger value for wider distribution spread)
kappa   = 0; %3 - L;          % from literature
beta    = 2;          % optimal in case distribution of x is Gaussian
lambda  = alpha^2*(L + kappa) - L;  % scaling factor for sigma points
c       = L + lambda;               % scaling factor for covariance
x_prior = zeros(L,1);
y_prior = zeros(M,1);
P_prior = zeros(L,L);
P_xy    = zeros(L,M);
R = sqrt(R);
Q = sqrt(Q);

% Weights
W_s(1) = lambda/c;
W_c(1) = lambda/c + (1 - alpha^2 + beta);
W_s (2:2*L+1) = 1/(2*c);  % for states
W_c (2:2*L+1) = 1/(2*c);  % for covariances

% Sigma points
X{1} = x;
xxx_t = x(:,ones(1,numel(x)));
X3 = [x xxx_t+S xxx_t-S]; 
for i = 2:2*L+1
    X{i} = X3(:,i);
end

%% Unscented transformation 
for i=1:2*L+1
    Xx{i}     = dfun(X{i},u); % sigma state propagation
    x_prior   = x_prior + W_s(i)*Xx{i};                                    % prior estimate
    Yy{i}     = mfun(Xx{i},u);                                             % sigma measuremnent
    y_prior   = y_prior + W_s(i)*Yy{i};                                    % measurement of prior
end
for i=1:2*L+1
    d_x(:,i)  = sqrt(abs(W_c(i)))*(Xx{i} - x_prior);
    d_y(:,i)  = sqrt(abs(W_c(i)))*(Yy{i} - y_prior);
    P_xy      = P_xy + W_c(i)*(Xx{i} - x_prior)*(Yy{i} - y_prior)';        % cross-covariance
end
%% S update
x_res = d_x;%*diag(sqrt(abs(W_c)));
y_res = d_y;%*diag(sqrt(abs(W_c)));
% Transpose of the matrix is used because qr returns an UT matrix;:
M_x = [x_res(:,2:2*L+1) Q]';
M_y = [y_res(:,2:2*L+1) R]';
[M_x_post] = givens_full(M_x); % qr([x_res(:,2:2*L+1) Q]',0);
[M_y_post] = givens_full(M_y); % qr([y_res(:,2:2*L+1) R]',0);
S_x = M_x_post(1:L,1:L);
S_y = M_y_post(1:M,1:M);
% Large alpha may result into the negative central weight
if W_c(1) < 0 
    S_x = cholupdate(S_x,x_res(:,1),'-');
    S_y = cholupdate(S_y,y_res(:,1),'-');
else
    S_x = cholupdate(S_x,x_res(:,1),'+');
    S_y = cholupdate(S_y,y_res(:,1),'+');
end
S_prior = S_x';
P_yy= S_y'*S_y;
%% Correction
K           = P_xy*pinv(P_yy);           % Kalman gain
x_posterior = x_prior + K*(y - y_prior); % State correction
U           = K*S_y;
for i=1:size(U,2)
     [~,r] = chol(S_x);
     if r ~= 0 || rank(S_x) ~= size(S_x,1)
         cc = max(max(abs(S_x)));
         S_x = S_x +eye(size(S_x))*cc*2;
     end
     S_x = cholupdate(S_x,U(:,i),'-');  % Cholesky update
end
S_posterior = S_x';
end


function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end

