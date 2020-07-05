%%
A = eye(3); 
B = [1; 0; 0];

A
dynfun = @(x,u) A*x + B*u;

x = [2; 2; 2];
u = 1;
result1 = dynfun(x, u)

%%
clear all; close all;
x = zeros(4,1);

L       = length(x);  % state vector dimension
% M       = length(y);  % measurement vector dimension
alpha   = 0.001;       % scaling constant - from literature (set to larger value for wider distribution spread)
kappa   = 3 - L;          % from literature
beta    = 2;          % optimal in case distribution of x is Gaussian
% lambda  = alpha^2*(L + kappa) - L;  % scaling factor for sigma points
% c       = L + lambda;   

w_0 = @(L,alpha)(alpha^2*(L + kappa) - L)/(L + alpha^2*(L + kappa) - L);
w_m = @(L,alpha)(1)/(2*(L + alpha^2*(L + kappa) - L));


a = [0.00001 0.0001 0.001 0.01 0.1 1];
ll = [3:10];

for j=1:length(ll)
    for i = 1:length(a)
        kappa = 0; %3 - ll(j);
        w(j,i) = w_0(ll(j),a(i));
        w_mm(j,i) = w_m(ll(j),a(i));
        summ(j,i) = w(j,i) + 2*ll(j)*w_mm(j,i);
    end
end

[a_x,ll_x] = meshgrid(ll,a);

figure;
surf(a_x,ll_x,w');
set(gca,'yscale','log')
set(gca,'zscale','log')
legend('w_0')

figure;
surf(a_x,ll_x,w_mm');
set(gca,'yscale','log')
set(gca,'zscale','log')
legend('w_m')

figure;
surf(a_x,ll_x,summ');
set(gca,'yscale','log')
% % set(gca,'zscale','log')
legend('sum of weights')

