%% exam1.m
% 1/9/2020
% ME327 - Exam 1
% Joshua Eckels
close all
clear all
clc
%% Problem 7
data = csvread('Exam1_data.csv');
x = data(:,1);
y = data(:,2);
n = size(x,1);
col1 = 1./x;
col2 = ones(n,1);
A = [col1 col2];
b = 1./y;
para = (A'*A)\(A'*b);
alpha = 1/para(2);
beta = para(1)*alpha;
y_est = alpha*(x./(beta+x)); % from equation 1

figure(1);
plot(x,y,'ok');
hold on
plot(x,y_est, '-k');
set(gcf,'color','w');
xlabel('X');
ylabel('Y (population)');
legend('Data', 'Regression model');

b_est = A*para;
e = b - b_est;
b_bar = mean(b);
total_sum = sum((b - b_bar).^2);
res_sum = sum(e.^2);
R2_lin = 1 - res_sum/total_sum;
fprintf("For alpha = %5.2f and beta = %5.2f, the R-squared value is %7.4f.\n",alpha,beta,R2_lin);