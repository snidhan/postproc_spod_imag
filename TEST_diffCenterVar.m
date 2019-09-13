% TEST  --  diffCenterVar
%
% Test diffCenterVar and demonstrate usage
%

clc; clear;

% Create a test function
a = 1.42;
b = 2.5;
c = 2.42;
xFun = @(t)( [a*cos(b*t + c); a^2*sin(c*t-b)] );
dxFun = @(t)( [-a*b*sin(b*t + c); a^2*c*cos(c*t-b)] );
ddxFun = @(t)( [-a*b*b*cos(b*t + c); -a^2*c*c*sin(c*t-b)] );

% Evaluate on a non-uniform grid
tBnd = [0,6];
nGrid = 100;
noiseScale = 0.2;
tGrid = linspace(tBnd(1), tBnd(2), nGrid);
tGrid = sort(tGrid + noiseScale*((tBnd(2)-tBnd(1))/nGrid)*randn(1,nGrid));
xGrid = xFun(tGrid);

% Derivatives!
[dxGrid, ddxGrid] = diffCenterVar(tGrid,xGrid);

% Truth model:
t = linspace(tBnd(1), tBnd(2), 200);
x = xFun(t);
dx = dxFun(t);
ddx = ddxFun(t);

% Plots:
figure(1523); clf;

h(1) = subplot(3,1,1); hold on;
plot(t,x)
plot(tGrid,xGrid,'o')
xlabel('time')
ylabel('position')
legend('truth','F.D.','Location','NorthEastOutside')
axis tight
title('Finite Differences, Non-Uniform Grid')

h(2) = subplot(3,1,2); hold on;
plot(t,dx)
plot(tGrid,dxGrid,'o')
xlabel('time')
ylabel('velocity')
legend('truth','F.D.','Location','NorthEastOutside')
axis tight


h(3) = subplot(3,1,3); hold on;
plot(t,ddx)
plot(tGrid,ddxGrid,'o')
xlabel('time')
ylabel('acceleration')
legend('truth','F.D.','Location','NorthEastOutside')
axis tight

linkaxes(h,'x');




