%% Eckels_Joshua_L10.m
% ME327 - Lab 10
% 2/21/2020
clear all
close all
clc

%% Lab 10
% parameters
k = 385;                % W/m-K
h = 40;                 % W/m^2-K
Tinf = 290;             % K
Tsurr = 250;            % K
eps = 0.8;              % emissivity
sig = 5.669*10^(-8);    % W/m^2-K^4
L = 0.40;               % m
D = 0.02;               % m
T0 = 400;               % K
Ac = (pi/4)*D^2;        % m^2
P = pi*D;               % m
param = [k,h,Tinf,Tsurr,eps,sig,L,D,T0,Ac,P];

R = @(vguess)residual(vguess,param);
dRdvg = @(vguess)deriv(vguess,R);
toggle = 1;
maxIter = 1000;
tol = 10^(-10);
vguess0 = 1000; % initial guess for shooting parameter dT/dx(x=0)

% Solve for shooting paramter with Newton-Raphson
[vfinal, n_iter] = newton_raphson(R,dRdvg,vguess0,tol,maxIter,toggle);

% Use ode45 to compute the temperature distribution using final IVs
absTol = 10^(-8);
relTol = 10^(-10);
options = odeset('AbsTol',absTol,'RelTol',relTol);
v0 = [T0;vfinal];    % Initial conditions
dvdx = @(x,v_vec)heat_deriv(x,v_vec,param);
[x_vec, v_mat] = ode45(dvdx, [0 L], v0,options);

% Display results
figure(1)
plot(x_vec,v_mat(:,1),'-r'); % plot Temperature distribution
xlabel('Distance along fin x [m]');
ylabel('Temperature [K]');
title('Temperature distribution along fin');

fprintf("Tip temperature of the fin is %6.4f K.\n",v_mat(end,1));

%% Functions

% Derivative function for heat transfer state system
% Parameters:
%   x     - location value to evaluate dvdx at
%   v_vec - state vector to evaluate dvdx at
% Outputs:
%   dvdx  - the output vector of derivatives of state vector
% Calling:
%   @(x,v_vec)heat_deriv(x, v_vec); params are passed in
function dvdx = heat_deriv(x, v_vec,param)
    % parameters passed in
    k = param(1);      % W/m-K
    h = param(2);      % W/m^2-K
    Tinf = param(3);   % K
    Tsurr = param(4);  % K
    eps = param(5);    % emissivity
    sig = param(6);    % W/m^2-K^4
    L = param(7);      % m
    D = param(8);      % m
    T0 = param(9);     % K
    Ac = param(10);    % m^2
    P = param(11);     % m
     
    n = size(v_vec,1); % v_vec comes in as column state vector
    dvdx = zeros(n,1); % dvdx will leave as column state vector
    dvdx(1) = v_vec(2);
    dvdx(2) = (h*P)/(k*Ac)*(v_vec(1)-Tinf) + (eps*sig*P)/(k*Ac)*(v_vec(1)^4 - Tsurr^4);
end

% Newton-Raphson
% Parameters:
%   R - residual function; R(x)
%   dRdx - derivative of residual function; dRdx(x)
%   xi - initial guess
%   tol - tolerance
%   maxIter - number of iterations maximum
%   toggle - 0; no print info. 1; print info
function [soln, n_iter] = newton_raphson(R, dRdx, xi, tol, maxIter, toggle)

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

% computes the residual of heat transfer system for vguess
% Parameters:
%   vguess - independent variable for shooting parameter
% Calling:
%   residual expected to be called as R(vguess); params are passed in
function res = residual(vguess,param)
    % parameters passed in
    k = param(1);      % W/m-K
    h = param(2);      % W/m^2-K
    Tinf = param(3);   % K
    Tsurr = param(4);  % K
    eps = param(5);    % emissivity
    sig = param(6);    % W/m^2-K^4
    L = param(7);      % m
    D = param(8);      % m
    T0 = param(9);     % K
    Ac = param(10);    % m^2
    P = param(11);     % m
    
    % Run ode45 to compute state vector at x=L using vguess
    absTol = 10^(-8);
    relTol = 10^(-10);
    options = odeset('AbsTol',absTol,'RelTol',relTol);
    v0 = [T0;vguess];    % Initial conditions
    dvdx = @(x,v_vec)heat_deriv(x,v_vec,param);
    [x_vec, v_mat] = ode45(dvdx, [0 L], v0,options);
    
    % Pull out end results from ode45 (at edge of fin)
    Tend = v_mat(end,1);
    dTdxend = v_mat(end,2);
    
    % Boundary Condition for dT/dx at x=L
    BC = (-h/k)*(Tend-Tinf) - (eps*sig/k)*(Tend^4 - Tsurr^4);
    
    % Compute residual
    res = dTdxend - BC;
end

% computes an approximate derivative of a residual function
% Parameters:
%   R - function handle to approximate derivative of
%   x - location to evaluate the derivative
% Calling:
% deriv expected to be called as dRdx(x); so R function handle must be
% passed in
function d = deriv(x,R)
    d = (R(1.01*x) - R(x))/(0.01*x);
end