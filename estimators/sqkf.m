function [x_new,x_prior,S_new,S_prior,Gain]  = sqkf(y,x,u,A,B,C,S,R,Q)
%Kalman Filter with square root Covariance update via
%triangularization
%Source: Simon, D. (2006): Optimal state estimation, Chapter 6, Page 169.
%And: Kaminski, P. (1971): Discrete Square Root Filtering:
%A Survey of Current Techniques, Page 731. 
Qsr = chol(Q,'lower');
Rsr = chol(R,'lower');
n = size(Qsr,2);
r = size(Rsr,2);
%Prediction
x_prior         = A*x + B*u;
prearray_prior  = [A*S, Qsr]';
postarray_prior = givens_full(prearray_prior);
S_prior         = postarray_prior(1:n,1:n)';
d_y      = y - C*x_prior;
prearray = [Rsr         , zeros(r,n);...
           (S_prior'*C'),  S_prior'];
%Correction
postarray   = givens_full(prearray); %or just %[~,postarray]= qr(prearray);
Kk          = postarray(1:r,r+1:(r+n))';
RR          = postarray(1:r,1:r);
S_new       = postarray(r+1:(n+r),r+1:(n+r))';
K           = (Kk/RR);
x_new       = x_prior + K * d_y;
Gain        = A - A*K*C;
% S_new = fcn_modGramSmith([S_new';Qsr']);
% S_new = S_new';
end

function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end

function W = fcn_modGramSmith(M)
%Source: Simon, D. (2006): Optimal state estimation, Chapter 6 ,Page 172.
n = size(M,2);
W = zeros(n);
for i=1:n
    sigma = sqrt(M(:,i)'*M(:,i));
    for j=1:n
        if j==i
            W(i,j)=sigma;
        elseif j==i+1 || j==n
            W(i,j)=(M(:,i)'*M(:,j))/sigma;
        else
            W(i,j)=0;
        end
    end
end
end
