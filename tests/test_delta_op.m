%%
my_init;
visFlag = 'On';
%%
dT = 0.01;
x = 1:dT:100;
y = sin(x);
y1 = cos(x);
y3 = -sin(x);
y5 = -cos(x);
%% Forward operator - order up to 3
for t=2:length(x)-1;
    y2(t-1) = delta_operator(1,y,t,dT,'forward');
end

for t=3:length(x)-2;
    y4(t-2) = delta_operator(2,y,t,dT,'forward');
end

for t=4:length(x)-3;
    y6(t-3) = delta_operator(3,y,t,dT,'forward');
end
%%
fig('Forward',visFlag);
subplot(3,1,1)
plot(y,'b'); hold on;
plot(y1,'k','LineWidth',2); hold on;
plot(y2,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 1');
subplot(3,1,2)
plot(y,'b'); hold on;
plot(y3,'k','LineWidth',2); hold on;
plot(y4,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 2');
subplot(3,1,3)
plot(y,'b'); hold on;
plot(y5,'k','LineWidth',2); hold on;
plot(y6,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 3');

%% Backward operator - order up to 3
for t=2:length(x)-1;
    yb2(t-1) = delta_operator(1,y,t,dT,'backward');
end

for t=3:length(x)-2;
    yb4(t-2) = delta_operator(2,y,t,dT,'backward');
end

for t=4:length(x)-3;
    yb6(t-3) = delta_operator(3,y,t,dT,'backward');
end
%%
fig('Backward',visFlag);
subplot(3,1,1)
plot(y,'b'); hold on;
plot(y1,'k','LineWidth',2); hold on;
plot(yb2,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 1');
subplot(3,1,2)
plot(y,'b'); hold on;
plot(y3,'k','LineWidth',2); hold on;
plot(yb4,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 2');
subplot(3,1,3)
plot(y,'b'); hold on;
plot(y5,'k','LineWidth',2); hold on;
plot(yb6,'r--','LineWidth',2); hold on;
legend('function','analytical','$\delta$-operator')
title('$\lambda$ = 3');

