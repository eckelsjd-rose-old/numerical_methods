%% lab4.m
% ME327 - Dr. Jones
% 1/10/20
% Joshua Eckels
clear all
close all
clc

%% Task 1
[sol, n] = func_MDnewton(@residual,@res_diff, [2 1]', 0.001, 10,1);
fprintf('Task 1-1: The final solution to 3 decimal places is: x1 = %5.3f, x2 = %5.3f.\n\n',sol(1),sol(2));
[sol, n] = func_MDnewton(@residual,@res_diff, [1 2]', 0.001, 10,1);
fprintf('Task 1-2: The final solution to 3 decimal places is: x1 = %5.3f, x2 = %5.3f.\n\n',sol(1),sol(2));

%% Task 2
% solve a general cubic function: (s-a)(s-b)(s-c)=0
% use initial guess [a,b,c]' = [1 3 5]' and an exact Jacobian
coeff = [-79.560 -1087.209 1046.133]'; % [A,B,C]'
[sol, n] = func_MDnewton(@(xi)residual2(xi,coeff),@res_diff2, [1 3 5]', 0.001, 100,1);
fprintf('Task 2: Using exact Jacobian\nthe three roots are: %7.3f, %7.3f, and %7.3f\n\n',sol(1),sol(2), sol(3));

%% Task 3
% solve a general cubic function: (s-a)(s-b)(s-c)=0
% use initial guess [a,b,c]' = [1 3 5]' and an approximate Jacobian
coeff = [-79.560 -1087.209 1046.133]'; % [A,B,C]'
[sol, n] = func_MDnewton(@(xi)residual2(xi,coeff),@(xi)approxJ(xi,coeff), [1 3 5]', 0.001, 100,1);
fprintf('Task 2: Using approximate\nthe three roots are: %7.3f, %7.3f, and %7.3f\n',sol(1),sol(2), sol(3));
