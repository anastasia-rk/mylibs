function[Q] = Qf(A,G,dT)
% m - scale structure, involves all info about scale
Q = zeros(size(A));
dt = dT/1000;
for t=0:dt:dT
    v = expm(A*(dT-t))*G;
    Q = Q + v*v'*dt;
end
end
