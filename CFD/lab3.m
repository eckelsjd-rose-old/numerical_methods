%% lab3.m
% ME427 - Intro to CFD
% Joshua Eckels
% 4/4/2020
clear all
close all
clc
%% Solve the Burgers equation for u(x), with nu = 0.1
% Start solution with N1 = 11 points. Then increment by 4 until we achieve
% an error < 0.001
fprintf("PART 1:\n");
norminf = 1;
tol = 0.001;    % desired convergence tolerance
N1 = 7;         % iteration starts at N1 + 4 = 11
nu = 0.1;
fprintf("Searching for solution to Burgers equation with diffusivity = %3.1f. . .\n",nu);
while (norminf > tol)
    % estimating error at N1; use N2 points for Richardson's extrapolation
    N1 = N1 + 4;
    N2 = 2*N1 - 1;
    % compute target discretization using N1 points
    [x1, urt1, dx1] = main_func(N1,nu);

    % compute discretization using N2 points
    [x2, urt2, dx2] = main_func(N2,nu);

    % approximate error in N1 discretization
    diff = zeros(1,size(x1,2));
    n = dx1/dx2; % amount dx was divided from N1 to N2
    for i = 1:size(urt1,1)
        diff(i) = urt2((i-1)*2+1) - urt1(i);
    end
    % Richardson's Extrapolation (2nd-order accuracy)
    error = (n^2/(n^2-1))*(diff);
    norminf = max(abs(error));
    fprintf("Error from %s at N=%d is estimated to be %6.4e m/s.\n","Finite difference",N1,norminf);
end

fprintf("Solution converged at N=%d. Error: %6.4e < %6.4e\n\n",N1,norminf,tol);
figure(1)
plot(x1,urt1','or');
xlabel('X coordinate [m]');
ylabel("Velocity 'u' [m/s]");
title(sprintf('Solution to 1D Burgers equation; N=%d, nu=%3.1f',N1,nu));
xlim([-1 1]);
ylim([-1 1]);
grid on
save ../lab4/fd.mat x1 urt1

%% Find lowest resolvable value of diffusivity
% fprintf("PART 2:\n");
% nu = 1;
% found_nu = false;
% N1 = 9;         % iteration starts at N1 + 2 = 11
% tmax = 7;       % time out after 7 seconds
% reduce = 0.1;   % reduce diffusivity by 10% each time
% while (~found_nu)
%     % Reduce diffusivity by 10% on each iteration
%     nu = (1-reduce)*nu;
%     
%     % Start solution with N1 from last iteration. Then increment by 2 and try to 
%     % achieve an error < 0.001. Quit when it takes longer than 7 seconds.
%     norminf = 1;
%     tol = 0.001;    % desired convergence tolerance
%     fprintf("Searching for solution to Burgers equation with diffusivity = %6.4e. . .\n",nu);
%     while (norminf > tol)
%         tic % start timing here
%         % estimating error at N1; use N2 points for Richardson's extrapolation
%         N1 = N1 + 2;
%         N2 = 2*N1 - 1;
%         % compute target discretization using N1 points
%         [x1, urt1, dx1] = main_func(N1,nu);
% 
%         % compute discretization using N2 points
%         [x2, urt2, dx2] = main_func(N2,nu);
% 
%         % approximate error in N1 discretization
%         diff = zeros(1,size(x1,2));
%         n = dx1/dx2; % amount dx was divided from N1 to N2
%         for i = 1:size(urt1,1)
%             diff(i) = urt2((i-1)*2+1) - urt1(i);
%         end
%         % Richardson's Extrapolation (2nd-order accuracy)
%         error = (n^2/(n^2-1))*(diff);
%         norminf = max(abs(error));
%         fprintf("Error from %s at N=%d is estimated to be %6.4e m/s.\n","Finite difference",N1,norminf);
%         
%         % break if solution took longer than tmax seconds to compute
%         % allow tolerance to be reached before breaking
%         tend = toc;
%         if (~found_nu && tend > tmax)
%             tout = tend;        % save timeout seconds
%             found_nu = true;
%         end
%     end
%     if (~found_nu)
%         fprintf("Solution converged at N=%d. Error: %6.4e < %6.4e\n\n",N1,norminf,tol);
%     end
% end
% fprintf("Timed out at %6.4e s.\n",nu,tout);
% fprintf("Lowest resolvable diffusivity: nu = %6.2e at N = %d\n",nu,N1);
% fprintf("Final error: %6.4e < %6.4e\n\n",norminf,tol);
% figure(2)
% plot(x1,urt1','or');
% xlabel('X coordinate [m]');
% ylabel("Velocity 'u' [m/s]");
% title(sprintf('Solution to 1D Burgers equation; N=%d, nu=%6.2e',N1,nu));
% xlim([-1 1]);
% ylim([-1 1]);
% grid on
% fprintf("END.\n");

%% Functions
% Estimate solution to 1D Burgers equation
% Parameters:   N   = number of mesh points
%               nu  = diffusivity
% Return:       x   = domain mesh points in x-direction [m]
%               urt = returned 'u' function values evaluated at x
%               dx  = uniform mesh spacing [m]
function [x, urt,dx] = main_func(N,nu)
    % Bounds of discrete problem
    xlow = -1;
    xupp = 1;
    dx = (xupp-xlow)/(N-1);
    x = xlow:dx:xupp;

    % Define residual and Jacobian functions
    param = [nu,dx];                % [diffusivity, mesh spacing]
    R = @(u)residual(u,x,param);
    J = @(u)jacobian(u,R);

    % Perform Newton-Raphson method on initial guess of u0
    maxIter = 1000;
    toggle = 0;
    tol = 0.0001;
    u0 = [1:N]';   % initial guess
    [urt, n_iter] = func_MDnewton(R,J,u0,tol,maxIter,toggle);
end

%% Functions
% computes a Residual vector for the given nonlinear DE
% Parameters:
%   u - vector of u values to compute residual at [u1;u2;u3;u4;u5...]
%   x - vector of discrete x locations at [x1,x2,x3,x4,x5...]
%   param - parameters needed to compute the residual: param = [nu,dx];
% Calling:
%   residual expected to be called as residual(u); x and param must be
%   provided
function res = residual(u,x,param)

    N = size(u,1);
    res = zeros(N,1);
    
    % Parameters passed in
    nu = param(1);  % diffusivity
    dx = param(2);  % mesh spacing [m]
    
    % residual of the middle discrete values
    for i = 2:N-1
        res(i) = (-nu/dx^2)*u(i+1) + (2*nu/dx^2)*u(i) + (-nu/dx^2)*u(i-1) + ... 
        (1/(4*dx))*u(i+1)^2 + (-1/(4*dx))*u(i-1)^2;
    end
    
    % boundary values
    res(1,1) = u(1) - 1;
    res(N,1) = u(N) + 1;
end

% computes an approximate Jacobian
% Parameters:
%   R - vector of residual functions; called as R(u)
%   u - location to evaluate jacobian at
% Calling:
%   jacobian expected to be called as jacobian(u); so R function handle must be
%   passed in
function ret = jacobian(u,R)

    N = size(u,1);
    ret = zeros(N);
    for i = 1:N
        u_new = u;
        u_new(i) = u(i)*1.01;
        ret(:,i) = (R(u_new) - R(u)) / (0.01*u(i));
    end
end

% Multi-dimensional Newton-Raphson method
% Parameters:
%   R - vector of residual functions; called as R(u)
%   J - approximate Jacobian function; called as J(u)
%   xi - initial location to evaluate at
%   tol - tolerance level to meet
%   maxIter - maximum number of iterations to perform
%   toggle - 0: hide print info, 1: display print info
function [xrt, n_iter] = func_MDnewton(R, J, xi, tol, maxIter, toggle)

    N = size(xi,1);         % size of input residual vector
    corr = -J(xi)\R(xi);    % correction factor
    xi_new = xi;
    n = 0; %                number of iterations
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
    if toggle == 1
         fprintf("Error from %s at N=%d is estimated to be %6.4e m/s.\n\n","Newton-Raphson",N,max(abs(corr)));
    end
end