function [Sx_p,S_p,Sx_pr,S_pr,K_p]  = sqrt_prior_to_prior(y,Sx,u,A,B,C,S,R,Q)
Qsr = chol(Q,'lower');
Rsr = chol(R,'lower');
n = size(Qsr,2);
r = size(Rsr,2);

M_pre = [Rsr C*S zeros(r,n);...
         zeros(r,n) A*S Qsr; ...
         zeros(r,n) eye(4) zeros(r,n);...
         -y'*Rsr' Sx' zeros(1,n)];   
%% unfinished
end

function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end