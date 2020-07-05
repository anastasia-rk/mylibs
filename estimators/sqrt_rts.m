function [x_sm,S_sm,Sx_sm_new,SS_new]  = sqrt_rts(Sx_sm,Sx_tt,Sx_tt1,SS_sm,SS_tt,S_tt1,Gain)
% square root RTS (stable version) from Park & Kaliath, 1995

% Gain   = (S_t)^(*/2)*(F_pt)'*(S_tt1)^(-*/2)
% S_tt1  = (P_tt1)^(1/2)
% Sx_sm  = (S_tT)^(-1)*x_tT - smoothed
% Sx_tt1 = (S_tt1)^(-1)*x_(tt1) -  prior from t times posterior from t
% Sx_tt  = (S_(t-1|t))^(-1)*x_tt - posterior
% SS_sm  = (S_tt1)^(-1)*(S_tT)*((S_tt1)^(-1)*(S_tT))'
% SS_tt  = (S_(t-1|t))^(-1)*(S_tt) - prior from t times posterior from t

Sx_sm_new = Sx_tt + Gain*(Sx_sm - Sx_tt1);
SS_new = SS_tt*SS_tt' - Gain*Gain' + Gain*SS_sm*Gain';

% recover the state and covariance
x_sm = S_tt1*Sx_sm_new;
P_sm = S_tt1*SS_new*S_tt1';
[LM,DM,E,pneg] = mchol(P_sm);
S_sm = LM*sqrt(DM);



