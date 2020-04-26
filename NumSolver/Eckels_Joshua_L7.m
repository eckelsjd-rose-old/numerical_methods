%% Eckels_Joshua_L7.m
% Joshua Eckels
% ME327 - Lab 7
% 31 Jan 2020
clear all
close all
clc
% read data
data = xlsread('acceleration-data.xls');
t = data(:,1);
a = data(:,2);
n = size(t,1);  % num data points
figure(1)
plot(t,a,'.k');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration v. time');

%% Rectangle Approximation
% find velocity at n/2 points with rectangle approx
vR = zeros(floor(n/2),1);
vtime = zeros(floor(n/2),1);
sum = 0;
for i = 1:floor(n/2)
    vtime(i,1) = t(2*i);
    dt = t(2*i+1,1) - t(2*i-1,1);
    sum = sum + a(2*i)*dt;
    vR(i,1) = sum;
end
% find position at n/4 points with rectangle approx
% since size(vR) is even, we should leave off last rectangle
% (we can't divide an even number evenly into groups of 3)
dR = zeros(floor(n/4)-1,1);
dtime = zeros(floor(n/4)-1,1);
sum = 0;
for i = 1:floor(n/4)-1
    dtime(i,1) = vtime(2*i);
    dt = vtime(2*i+1,1) - vtime(2*i-1,1);
    sum = sum + vR(2*i)*dt;
    dR(i,1) = sum;
end
figure(2)
subplot(211);
plot(vtime,vR,'.k');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Rectangle approx. of velocity');
subplot(212);
plot(dtime,dR,'.k');
xlabel('Time [s]');
ylabel('Position [m]');
title('Rectangle approx. of position');

%% Trapezoid Approximation
% find velocity at n-1 points with trapezoid approx
vR = zeros(n-1,1);
vtime = zeros(n-1,1);
sum = 0;
for i = 1:n-1
    vtime(i,1) = t(i);
    dt = t(i+1,1) - t(i,1);
    sum = sum + ((a(i)+a(i+1))/2)*dt;
    vR(i,1) = sum;
end
Tapprox1200 = sum;
% find position at n-2 points with trapezoid approx
dR = zeros(n-2,1);
dtime = zeros(n-2,1);
sum = 0;
for i = 1:n-2
    dtime(i,1) = vtime(i);
    dt = vtime(i+1,1) - vtime(i,1);
    sum = sum + ((vR(i)+vR(i+1))/2)*dt;
    dR(i,1) = sum;
end
figure(3)
subplot(211);
plot(vtime,vR,'.k');
xlabel('Time [s]');
ylabel('Instantaneous Velocity [m/s]');
title('Trapezoid approx. of velocity');
subplot(212);
plot(dtime,dR,'.k');
xlabel('Time [s]');
ylabel('Instantaneous Position [m]');
title('Trapezoid approx. of position');

%% Task 5
% Redo the Trapezoid approximation of velocity when N=600
% Then compute the estimated error when N=1200
vR = zeros(floor(n/2),1);
vtime = zeros(floor(n/2),1);
sum = 0;
for i = 1:floor(n/2)
    vtime(i,1) = t(2*i);
    dt = t(2*i+1,1) - t(2*i-1,1);
    sum = sum + ((a(2*i-1)+a(2*i+1))/2)*dt;
    vR(i,1) = sum;
end
Tapprox600 = sum;
error1200 = (Tapprox1200 - Tapprox600)/3;
fprintf(['Trapezoid method: the error estimate for the vehicle speed at',...
' the last time record is %9.2e m/s\n'],error1200);

%% Task 6
% Use regression with cubic function to get surrogate model
A = [ones(size(dtime,1),1) dtime./100 (dtime./100).^2 (dtime./100).^3];
b = dR;
beta = (A'*A)\(A'*b);
d_est = A*beta;
figure(4);
plot(dtime,d_est,'-r',dtime,dR,'.k');
legend('Surrogate model','Data');
xlabel('Time [s]');
ylabel('Instantaneous Position [m]');
title('Trapezoid approx. of position');
fprintf(['Regression:\nBeta0=%5.2f Beta1=%5.2f Beta2=%5.2f Beta3=%5.2f\n'...
    ''],beta(1),beta(2),beta(3),beta(4));

%% Task 7
tol = 0.02;
toggle = 0;
maxIter = 1500;
% beta is declared above from linear regression

dist = [5000 10000 15000]'; % distances
t_init = [500 800 1000]';   % initial guesses
tresult = zeros(3,1);       % result times
iter_result = zeros(3,1);   % result iteration count
for i = 1:3
    t0 = t_init(i);
    % arr_time called like arr_time(distance), expects all other vars 
    arr_time = @(d)main_func(d,tol,toggle,maxIter,t0,beta);
    [tresult(i),iter_result(i)] = arr_time(dist(i));  
    fprintf("It takes %7.2f seconds for the vehicle to travel %2i km.\n",tresult(i),dist(i));
end
fprintf("END\n");

%% Functions %%
% Computes arrival times given a distance d in km
% Gets passed all arguments needed by residual,deriv, and func_newtwon
% functions, and calls them as needed
function [time, iter] = main_func(d,tol,toggle,maxIter,t0,beta)
    % rh will be called like rh(t), and expects d and beta to exist
    rh = @(t)residual(t,d,beta);        % residual function handle
    % dh will be called like dh(t), and expects rh(t) to exist
    dh = @(t)deriv(t,rh);               % derivative function handle
    [time, iter] = func_newton(rh,dh,t0,tol,maxIter,toggle);
end

% Solves a residual function using Newton-Raphson method
% Parameters:
%   R - residual function
%   dRdx - derivative of residual function
%   xi - initial guess
%   tol - error tolerance
%   maxIter- maximum number of iterations to run before stopping
%   toggle - 1: show print info, 0: hide print info
function [soln, n_iter] = func_newton(R, dRdx, xi, tol, maxIter, toggle)
corr = abs(R(xi)/dRdx(xi));
xi_new = xi;
n = 0; % number of iterations
if (toggle == 1)
    fprintf('%5s %7s %9s %9s %9s\n','Count','xi','R(xi)','dRdx(xi)','corr');
end

while (abs(corr) > tol)
    n = n + 1;
    Rxi = R(xi);
    dRdxi = dRdx(xi);
    corr = -Rxi/dRdxi;
    xi_new = xi + corr;
    if (toggle == 1)
        fprintf('%5.0f %7.3f %9.3f %9.3f %9.3f\n', n, xi, Rxi, dRdxi, corr);
    end
    xi = xi_new;
    if (n >= maxIter)
        break;
    end
end
% found solution (or exceeded maxIter)
soln = xi_new;
n_iter = n;
end

% computes the residual using cubic surrogate of a model
% Parameters:
%   t - independent variable for time
%   beta - cubic model coefficients in an array
%   d - desired function value
% Calling:
%   residual expected to be called as R(t); so d and beta must be provided
function res = residual(t,d,beta)
    B0 = beta(1);
    B1 = beta(2);
    B2 = beta(3);
    B3 = beta(4);
    res = B0 + B1*(t/100) + B2*(t/100)^2 + B3*(t/100)^3 - d;
end

% computes an approximate derivative
% Parameters:
%   R - function handle to approximate derivative of
%   t - location to evaluate the derivative (time)
% Calling:
% deriv expected to be called as dRdx(t); so R function handle must be
% passed in
function d = deriv(t,R)
    d = (R(1.01*t) - R(t))/0.01;
end