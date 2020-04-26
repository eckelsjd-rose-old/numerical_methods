%% P9_Eckels.m
% Joshua Eckels
% ME327 - Final exam
% 2/26/2020
clear all
close all
clc
%% Problem 9
v0 = [1;2;3];       % Initial conditions
dvdx = @(x,v_vec)state_deriv(x,v_vec);
xf = 10;            % simulate to x=10
dx = 0.001;         % step of 0.1 m
N = xf/dx;          % number of steps
[x_vec, v_mat] = rk2(dvdx,v0,dx,N);

figure(2)
grid on
plot(x_vec,v_mat(:,1),'.k');
xlabel('Distance x [m]');
ylabel('Function value f [-]');
title('Problem 9: RK2 approximation');

fprintf("Problem 9:\n");
fprintf("At x=%d, derivative dfdx=%6.4f\n",x_vec(end),v_mat(end,2));

% approximate the error with Richardson's Extrapolation
% result from N=10,000 saved above, repeat with N=20,000
dx2 = dx/2;         % divide step in half
N2 = xf/dx2;        % number of time steps
[x_vec2, v_mat2] = rk2(dvdx,v0,dx2,N2);

% find error in dfdx from N=10,000 to N=20,000 
xdiff = zeros(1,size(x_vec,1));
n = dx/dx2; % amount dx was divided by from N=10,000 to N=20,000
for i = 1:size(v_mat,1)
    xdiff(i) = v_mat2(((i-1)*2+1),2) - v_mat(i,2);
end
x_error = (n^2/(n^2-1))*(xdiff); % second order accurate for rk2
norminf = max(abs(x_error));

fprintf("Error in function value f at x=%d m is estimated to be %6.8e.\n", x_vec(end),norminf);

%% Functions
% Runge-Kutta approximation with 2 slopes
% Parameters:
%   dv_fcn - function handle which evaluates derivative of state vector
%   v0_fcn - vector with initial value of each state variable [nx1]
%   dx     - time step 
%   N      - number of steps
% Outputs:
%   x_vec  - vector of x values [ (N+1) x 1 ]
%   v_mat  - matrix with instantaneous value of each state variable 
%            [ (N+1) x n ]
% Calling:
%   dv_fcn expected to be called as dv_fcn(x, v_vec), where v_vec is the
%   current state vector
function [x_vec, v_mat] = rk2(dv_fcn, v0_vec, dx, N)
    n = size(v0_vec,1);     % n = number of state variables
    x_vec = (0:dx:dx*N)';   % x column vector
    v_mat = zeros(N+1,n);   % columns for each state var; rows for time steps
    v_mat(1,:) = v0_vec';   % row=0; x=0; v = v0 (ICs)
    for i = 1:N
        k1 = dv_fcn(x_vec(i),v_mat(i,:)'); % d{v}/dt 1 (col vector)
        ystar1 = v_mat(i,:)' + (dx/2)*k1;
        
        k2 = dv_fcn(x_vec(i) + dx/2, ystar1); % d{v}/dt 2
        v_mat(i+1,:) = v_mat(i,:) + dx*k2';
    end
end

% Derivative function for the differential equation in Problem 9
% Parameters:
%   x     - x value to evaluate dvdx at
%   v_vec - state vector to evaluate dvdx at
% Outputs:
%   dvdx  - the output vector of derivatives of state vector
% Calling:
%   @(x,v_vec)pend_deriv(x, v_vec)
function dvdx = state_deriv(x, v_vec)
     n = size(v_vec,1); % v_vec comes in as column state vector
     dvdx = zeros(n,1); % dvdt will leave as column state vector
     dvdx(1) = v_vec(2);
     dvdx(2) = v_vec(3);
     dvdx(3) = -(1/2)*v_vec(1)^2 * v_vec(2);
end