%% Eckels_Joshua_L8.m
% ME327 - Lab 8
% Joshua Eckels
% 2/7/2020
clear all
close all
clc

%% Task 3
% compute disretization approximation when N=11
[x11, yrt11, dx11] = main_func(11);

% compute disretization approximation when N=21
[x21, yrt21, dx21] = main_func(21);

% find error from N=11 to N=21
diff = zeros(1,size(x11,2));
n = dx11/dx21; % amount dx was divided by from N=11 to N=21
for i = 1:size(yrt11,1)
    diff(i) = yrt21((i-1)*2+1) - yrt11(i);
end
error = (n^2/(n^2-1))*(diff);
norminf = max(abs(error));
fprintf("Error from %s at N=%d is estimated to be %6.4e m.\n","Finite difference",11,norminf);

figure(1)
plot(x11,yrt11'*100,'or');
xlabel('X coordinate [m]');
ylabel('Beam deflection [cm]');
title('Beam deflection along half-length');
xlim([0 2]);
ylim([-50 10]);
grid on
fprintf("END.\n");

%% Functions
% Estimate solution to beam deflection problem using N points at initial
% guess of y0
function [x,yrt,dx] = main_func(N)
    % Parameters
    E = 10*10^9;        % Elastic modulus [Pa]
    b = 0.20;           % Beam width [m]
    h = 0.10;           % Beam height [m]
    I = (1/12)*b*h^3;   % Moment of inertia [m^4]
    L = 2;              % Beam half-length [m]
    P = -50*10^3;       % Applied load [N]

    % Bounds of discrete problem
    xlow = 0;
    xupp = L;
    dx = (xupp-xlow)/(N-1);
    x = xlow:dx:xupp;

    % Define residual and Jacobian functions
    param = [E,I,P,dx];
    R = @(y)residual(y,x,param);
    J = @(y)jacobian(y,R);

    % Perform Newton-Raphson method on initial guess of y0
    maxIter = 1000;
    toggle = 0;
    tol = 0.0001;
    y0 = [1:N]';   % initial guess
    [yrt, n_iter] = func_MDnewton(R,J,y0,tol,maxIter,toggle);
    fprintf("The maximum deflection is %5.2f cm at an axial location of %5.2f m for N=%d.\n", yrt(N)*100,x(N),N);
end

% computes a Residual vector for the given nonlinear DE
% Parameters:
%   y - vector of y values to compute residual at [y1;y2;y3;y4;y5...]
%   x - vector of discrete x locations at [x1,x2,x3,x4,x5...]
%   param - parameters needed to compute the residual: param = [E,I,P,dx];
% Calling:
% residual expected to be called as residual(y); x and param must be
% provided
function res = residual(y,x,param)

    N = size(y,1);
    res = zeros(N,1);
    
    % Parameters passed in
    E = param(1);   % Modulus of Elasticity [Pa]
    I = param(2);   % Moment of inertia [m^4]
    P = param(3);   % Applied load [N]
    dx = param(4);  % delta x [m]
    
    % residual of the middle discrete values
    for i = 2:N-1
        res(i) = (y(i+1)-2*y(i)+y(i-1))/(dx^2) + ((P*x(i))/(2*E*I))*(1+((y(i+1)-y(i-1))/(2*dx))^2)^(3/2);
    end
    
    % boundary values
    res(1,1) = y(1);
    res(N,1) = ((1/2)*y(N-2)-2*y(N-1)+(3/2)*y(N))/dx;
end

% computes an approximate Jacobian
% Parameters:
%   R - vector of residual functions; called as R(y)
%   y - location to evaluate jacobian at
% Calling:
% jacobian expected to be called as jacobian(y); so R function handle must be
% passed in
function ret = jacobian(y,R)

    N = size(y,1);
    ret = zeros(N);
    for i = 1:N
        y_new = y;
        y_new(i) = y(i)*1.01;
        ret(:,i) = (R(y_new) - R(y)) / (0.01*y(i));
    end
end

% Multi-dimensional Newton-Raphson method
% Parameters:
%   R - vector of residual functions; called as R(y)
%   J - approximate Jacobian function; called as J(y)
%   xi - initial location to evaluate at
%   tol - tolerance level to meet
%   maxIter - maximum number of iterations to perform
%   toggle - 0: hide print info, 1: display print info
function [xrt, n_iter] = func_MDnewton(R, J, xi, tol, maxIter, toggle)

    N = size(xi,1); % size of input residual vector
    corr = -J(xi)\R(xi);
    xi_new = xi;
    n = 0; % number of iterations
    if (toggle == 1)
        fprintf('%9s %12s   %1s\n','Iteration','Max Error','Current root approximation');
    end

    while (max(abs(corr)) > tol)
        n = n + 1;
        Rxi = R(xi);
        Jxi = J(xi);
        corr = -Jxi\Rxi;
        xi_new = xi + corr;
        if (toggle == 1)
            g = sprintf('%9.3f', xi);
            fprintf('%6i    %12.2e %s\n', n, max(abs(corr)), g);
        end
        xi = xi_new;
        if (n >= maxIter)
            break;
        end
    end
    % found solution (or exceeded maxIter)
    xrt = xi_new;
    n_iter = n;
    fprintf("Error from %s at N=%d is estimated to be %6.4e m.\n","Newton-Raphson",N,max(abs(corr)));
end