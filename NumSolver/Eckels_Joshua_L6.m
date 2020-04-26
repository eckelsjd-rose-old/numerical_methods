%% lab6.m
% ME327
% Joshua Eckels
% 1/24/2020
clear all
close all
clc

% parameters
Lp = 43.2;          % cm
dw = 5;             % cm
ds = 2.5;           % cm

% constraints on Lwcg
lb = (ds+dw)/2;     % lower bound [cm]
ub = Lp-(ds+dw)/2;  % upper bound [cm]

% experimental data
x = [4 10 20 30 39 10.52]';
tswing = [0.242 0.230 0.244 0.276 0.304 0.23]';
% original
% x = [4 10 20 30 39]';
% tswing = [0.242 0.230 0.244 0.276 0.304]';

figure(1)
plot(x,tswing,'ok');
set(gcf,'color','w');
xlabel('Location of mass Lwcg [cm]');
ylabel('Residual pendulum swing time [s]');

% linear regression
A = [ones(size(x)) x x.^2 x.^3];
b = tswing;
beta = (A'*A)\(A'*b);
t_est = A*beta;
hold on
plot(x,t_est,'-r');
legend('Experimental','Regression');
fprintf("Regression:\nBeta0=%4.4f Beta1=%4.4f Beta2=%4.4f Beta3=%4.4f\n",beta(1),beta(2),beta(3),beta(4));

% plot fit
xVec = linspace(lb,ub,101)';
t_fit = surrogate_model(xVec,beta);

figure(2)
plot(xVec, t_fit,'-r');
hold on
plot(x,tswing,'ok');
legend('Surrogate','Experimental Data');

% optimize surrogate
[x_opt,fval] = fmincon(@(x)surrogate_model(x,beta),12,[],[],lb,ub);
fprintf("The optimal adjustable distance is %5.2f cm and \nthe minimum swing time is %5.3f seconds.\n",x_opt,fval);
% blackbox_simulator(x_opt);

function [y] = surrogate_model(x,beta)

B0 = beta(1);
B1 = beta(2);
B2 = beta(3);
B3 = beta(4);

y = B0 + B1*x + B2*x.^2 + B3*x.^3;
end