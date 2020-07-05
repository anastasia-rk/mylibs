function [x_new,S_posterior,Gain]  = sqrt_meas(y,x_prior,A,C,S_prior,S_pr_x,R)
Rsr = chol(R,'lower');
n = length(x_prior);
r = size(Rsr,2);
%Prediction
d_y         = y - C*x_prior;
prearray    = [        Rsr, zeros(r,n), zeros(n,n), -Rsr*y;...
               S_prior'*C',   S_prior',     eye(n), S_pr_x];
%Correction
postarray   = givens_full(prearray); %or just %[~,postarray]= qr(prearray);
Kk          = postarray(1:r,r+1:(r+n))';
RR          = postarray(1:r,1:r);
S_posterior = postarray(r+1:(n+r),r+1:(n+r))';
S_pr_S_post = postarray();
S_post_x_post = postarray();
K           = (Kk/RR);
Gain        = A - A*K*C;
end

function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end