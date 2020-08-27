 my_init;
 %%
% Time and scale;
l = 4; % state vector size
r = 2; % measurement vector size;
M = 6; % finest scale order;
dt = 0.001; % sstep on finest scale;
t0 = 0; % start time on finest scale;
tf = (20-M)*dt*2^M - dt; % end time on finest scale;

% Defining scales
for i=1:M+1
    m = i-1;
    scale(i).m = m;
    scale(i).dt = dt*2^(M-m);
    scale(i).dt_up = scale(i).dt/2;
    scale(i).dt_down = scale(i).dt/4;
    scale(i).N = M-m;
end

for i=M+1:-1:1
   % set starting and ending time
    sm = 0;
    if (i < M+1)
        for j = i:M
            sm = sm + scale(j+1).dt_up;
        end
    end
    scale(i).t0 = t0 + sm;
    scale(i).tf = tf - sm;
    
    % Placeholders for states and measurements
    scale(i).t = [scale(i).t0:scale(i).dt:scale(i).tf];
    if isempty(scale(i).t)
        scale(i).t = scale(i).t0;
    end
    X{i} = zeros(l,length(scale(i).t));
    Y{i} = zeros(r,length(scale(i).t));
    U{i} = zeros(1,length(scale(i).t));
    Scale_plot{i} = (M+1 - scale(i).m)*ones(1,length(scale(i).t));
end

%% Check scales and times
figure;
for i=1:M+1
    stem(scale(i).t,Scale_plot{i},'LineWidth',2); hold on;
end
grid on;
xlim([0,tf]);
xlabel('time,s'); ylabel('scale');
for i=1:M+1
    names(i) = {['m = ' num2str(M+1 - i)]};
end
set(gca,'ytick',[1:M+1],'yticklabel',names);
%% Integrate the system
l = 4; % state vector size
eps = 0.01; % mass ratio
m1 = 1;    % first mass
kk = 98.1;  % first spring tension
k1 = kk;

m2 = m1*eps;
k2 = kk;

% Noises
I = eye(l);
R_wsd = 10^5; 
Q     = 10^(-6)*I;
R = zeros(r);
dti = 0.0001;
mu  = 1;

% System matrices
A = [  0      1    0         0;...
    -k1/m1    0    k1/m1     0;...
       0      0    0         1;...
     k1/m2    0  -(k1+k2)/m2 0];
 
B = [0; 0; 0; k2/m2]; 
C = [1 0 0 0; 0 0 1 0];
D = 0;
G = [0; 0; 0; 1];

% To generate u
A_u = -mu;
B_u = 1;
C_u = 1;

% State space
sys = ss(A,B,C,D);
sys_noise = ss(A,G,C,D);
sysu = ss(A_u, B_u, C_u, 0);
sysud = c2d(sysu,dti);
sysd = c2d(sys,dti);

% Solve ssm
time1 = [t0:dti:tf]; 
w = normrnd(0,1,[1,length(time1)]);
u = lsim(sysud,w,time1,0);
[y,time,x] = lsim(sysd,u,time1,zeros(l,1));
y = y';
x = x';
u = u';

%% Smooth the origial system - both kf and rts are working
    X1(:,1) = zeros(l,1);
    P1(:,:,1) = eye(l);
    X_pr1(:,1) = zeros(l,1);
    Qq = Qf(sysd.a,G,dti);
    for t = 2:length(time1)
        [X1(:,t),P1(:,:,t)] = kf(y(:,t),X1(:,t-1),u(:,t-1),P1(:,:,t-1),sysd.a,sysd.b,C,Qq,R);
    end  
    T = t;
    X_sm1(:,T) = X1(:,T);
    P_sm1(:,:,T) = P1(:,:,T);
    for t = T-1:-1:1
        [X_sm1(:,t),P_sm1(:,:,t)] = rts(X_sm1(:,t+1),X1(:,t),P_sm1(:,:,t+1),P1(:,:,t),u(:,t),sysd.a,sysd.b,Qq);
    end  
    
    figure;
    plot(time,y(1,:),'k'); hold on;
    plot(time1,X_sm1(1,:),'k'); hold on;
    clear sysd sysd_noise
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Assign measurements to certain scales
% Sample the signal
K_r = 10^(4);
K_r1 = 10^(4);%10^10;
K_r2 = 10^(4);%10^10;
for i=1:M+1
    j = 1;
    for t = scale(i).t
        match = find(abs(time - t) <= 10^(-6));
        X{i}(:,j) = x(:,match);
        U{i}(:,j) = u(:,match);
        j = j + 1;
    end 
    % Discrete system matrices
    sysd = c2d(sys,scale(i).dt);
    sysd_noise = c2d(sys_noise,scale(i).dt);
    scale(i).A = sysd.a;
    scale(i).B = sysd.b;
    scale(i).C = zeros(size(C));
    scale(i).G = sysd_noise.b;
    % G actually should be eye(l)
    scale(i).Q = Qf(A,G,scale(i).dt);
    scale(i).R = K_r*eye(size(C,1));
    % Fine-to-coarse matrices
    sysd_f = c2d(sys,scale(i).dt_up);
    sysd_noise_f = c2d(sys_noise,scale(i).dt_up);
    scale(i).A_up = sysd_f.a;
    scale(i).B_up = sysd_f.b;
    scale(i).G_up = sysd_noise_f.b;
    scale(i).Q_up = Qf(A,G,scale(i).dt_up);
    % Coarse-to-fine matrices
    sysd_d = c2d(sys,scale(i).dt_down);
    sysd_noise_d = c2d(sys_noise,scale(i).dt_down);
    scale(i).A_down = sysd_d.a;
    scale(i).B_down = sysd_d.b;
    scale(i).G_down = sysd_noise_d.b;
    scale(i).Q_down = Qf(A,G,scale(i).dt_up);
    clear sysd sysd_noise sysd_f sysd_noise_f sysd_d sysd_noise_d
end

% Only measurements at certain scales are available
scls = [4];
for i=scls
 j = 1;
    for t = scale(i).t
        match = find(abs(time - t) <= 10^(-6));
        Y{i}(:,j) = y(:,match);
        j = j + 1;
        scale(i).C = C;
        scale(i).R = R;
    end 
end
% Plot samples;
fig('States','On');
subplot(2,1,1)
plot(time,y(1,:),'k'); hold on;
for i=1:M+1
        plot(scale(i).t,Y{i}(1,:),'-'); hold on;
        names1(i) = {['m = ' num2str(i-1)]};
end
legend(['noiseless measurement',names1]);
xlabel('time,s'); ylabel('$x_1$');
subplot(2,1,2)
plot(time,y(2,:),'k'); hold on;
for i=1:M+1
        plot(scale(i).t,Y{i}(2,:),'-'); hold on;
end
legend(['noiseless measurement',names1]);
xlabel('time,s'); ylabel('$x_2$');

fig('States','On');
plot(time,y(1,:),'k'); hold on;
for i=1:M+1
        plot(scale(i).t,Y{i}(1,:),'*','Linewidth',1); hold on;
end
legend(['noiseless measurement',names1]);
xlabel('time,s'); ylabel('$x_1$');

figure; 
plot(time,y(2,:),'k'); hold on;
for i=1:M+1
        plot(scale(i).t,Y{i}(2,:),'*','Linewidth',1); hold on;
end
legend(['noiseless measurement',names1]);
xlabel('time,s'); ylabel('$x_2$');



