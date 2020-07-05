function [x_prior,S_prior]  = sqrt_predict(y,x,u,A,B,S,Q)
%Kalman Filter with square root Covariance update via
%triangularization
%Source: Simon, D. (2006): Optimal state estimation, Chapter 6, Page 169.
%And: Kaminski, P. (1971): Discrete Square Root Filtering:
%A Survey of Current Techniques, Page 731. 
Qsr = chol(Q,'lower');
n = size(Qsr,2);
%Prediction
x_prior         = A*x + B*u;
prearray_prior  = [A*S, Qsr]';
postarray_prior = givens_full(prearray_prior);
S_prior         = postarray_prior(1:n,1:n)';
end

function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end