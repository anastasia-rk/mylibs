function [x_sm,S_sm,Sx_sm]  = sqrt_rts_fast(x_t1T,S_t1T,x_tt,S_tt,u,dfun,A,Q)
% square root RTS for invertable transition matrix
n          = length(x_tt);
x_t1t      = dfun(x_tt,u_p);
P_t1t(:,:) = A*S_tt*S_tt'*A' + Q;
Q_tilde    = Q - Q*pinv(P_t1t)*Q;
Q_sr       = chol(Q_tilde,'lower');
Gain       = S_tt*S_tt'*A'*pinv(P_t1t);
M_pre      = [Gain*S_t1T, pinv(A)*Q_sr]';
M_post     = givens_full(M_pre);
S_sm       = M_post(1:n,1:n);

Sx_sm_new = Sx_tt + Gain*(Sx_sm - Sx_tt1);
SS = S_t_S_tt*S_t_S_tt' - Gain*Gain' + Gain*SS_sm*SS_sm'*Gain';

% recover the state and covariance
x_sm = S_pr*Sx_sm_new;

end

function [M] =  givens_full(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = givens(M,i,j);
        end
    end
end


