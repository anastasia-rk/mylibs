clc; clear; close all;
t0 = 0; % start time on finest scale;

% Description of a double-mass double-spring system
eps = 0.01; % mass ratio
m1 = 1;    % first mass
kk = 98.1;  % first spring tension
k1 = kk;

m2 = m1*eps;
k2 = kk;

% Noise
dti = 0.01;
mu  = 1;
I = eye(4);
R_wsd = 10^(5); 
Q     = 10^(-4)*I;

% System matrices
A =   [  0      1    0         0;...
      -k1/m1    0    k1/m1     0;...
         0      0    0         1;...
       k1/m2    0  -(k1+k2)/m2 0];
 
B = [0; 0; 0; k2/m2];
C = eye(4);
D = 0;
G = sqrt(R_wsd)*[0; 0; 0; 1];

% To generate u
A_u = -mu;
B_u = 1;
C_u = 1;

% State space
sys = ss(A,B,C,D);
sysu = ss(A_u, B_u, C_u, 0);
sysud = c2d(sysu,dti);
sysd = c2d(sys,dti);

% Solve ssm
time = [t0:dti:10]; 
w = normrnd(0,1,[1,length(time)]);
u = lsim(sysud,w,time,0);
[y,t,x] = lsim(sysd,u,time,zeros(4,1));

%%
figure; set(gcf,'color','w');
subplot(3,1,1);
plot(t,x(:,1));
xlabel('$t$, s');ylabel('$x_1$');
subplot(3,1,2);
plot(t,x(:,3));
xlabel('$t$, s');ylabel('$x_2$');
subplot(3,1,3);
plot(t,u);
xlabel('$t$, s');ylabel('$u$');
ChangeInterpreter(gcf,'Latex');
