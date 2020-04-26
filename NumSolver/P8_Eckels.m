%% P8_Eckels.m
% Joshua Eckels
% ME327 - Final exam
% 2/26/2020
clear all
close all
clc
%% Problem 8
tlow = 0;
tf = 10;
dt = 0.001;
N = (tf-tlow)/dt + 1;
A = zeros(N);
b = zeros(N,1);
t = tlow:dt:tf;

% boundary conditions
% x(t=0) = 1; position N-1
b(N-1) = 1;
A(N-1,1) = 1;

% dxdt(t=0) = 0; position N
b(N) = 0;
A(N,1) = -1/dt;
A(N,2) = 1/dt;

% position 1:N-2
for i = 1:N-2
    A(i,i) = (2/dt^2) - 0.5/dt + 5;
    A(i,i+1) = (-4/dt^2) + 0.5/dt;
    A(i,i+2) = 2/dt^2;
    b(i) = exp(-t(i));
end
x = A\b;
figure(1)
plot(t,x,'.k');
xlabel('Time t [s]');
ylabel('Distance x [m]');
str = sprintf("Problem 8: Discretization using N=%d",N);
title(str);
fprintf("The value of x at tf=%d is %6.4f\n",tf,x(end));