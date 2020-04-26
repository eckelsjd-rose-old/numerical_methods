%% lab5.m
% Joshua Eckels
% ME427 - Intro to CFD
% Lab 5
% 4/20/20
clc
close all
clear all

%% Problem 1: Transverse mesh (y-direction)
% Tuning parameters for fGS = -A*exp(alpha*s)
N = 50;
alpha = 3;
ymin = 0;
ymax = -20;
A = ((ymin-ymax)*alpha)/(exp(alpha)-1);

% computational domain
ds = 1/(N-1);            % domain spacing
s = 0:ds:1;                 % domain range
fGS = -A*exp(alpha*s);       % grid-spacing function fGS

% physical domain
y = ymin-(A/alpha)*exp(alpha*s)+(A/alpha);
dy = y(2:length(y)) - y(1:length(y)-1);
gsr = (dy(2:length(dy)) - dy(1:length(dy)-1))./dy(1:length(dy)-1);
fprintf("fGS mesh results for N = %d:\n",length(y));
fprintf("dy max   = %4.3f\n",max(dy));
fprintf("dy min   = %4.3f\n",min(dy));
fprintf("max/min  = %4.3f\n",max(dy)/min(dy));
fprintf("max gsr  = %4.3f%%\n",max(gsr)*100);

% plot results
figure() % plot fGS(s)
subplot(2,1,1) 
plot(s,fGS,'-k');
str = sprintf("fGS Y mesh; N = %d",N);
xlabel('Computational domain [s]');
ylabel('Grid-spacing [fGS]');
title(str);
subplot(2,1,2) 
plot(s,y,'or');
xlabel('Computational domain [s]');
ylabel('Physical domain [y]');
title(str);
visualize_xy_mesh(0,1,ymin,ymax,[0,1],y);

% Transverse mesh (y-direction) using fGSR for grid-stretching ratio
A = -1/8;                    % fGSR tuning parameters
b = 2;                       %
c = 1/2;                     %
param = [ymin,ymax,N,A,b,c];

R = @(vguess)residual(vguess,param);    % Residual
dRdvg = @(vguess)deriv(vguess,R);       % Residual derivative
toggle = 0;
maxIter = 1000;
tol = 10^(-6);
vguess0 = 0.01; % initial guess for shooting parameter dy/ds(s=0)

% Solve for shooting paramter with Newton-Raphson
[vfinal, n_iter] = newton_raphson(R,dRdvg,vguess0,tol,maxIter,toggle);

% Use rk4 to compute the y-direction mesh from s=0 to s=1
v0 = [ymin;vfinal];                     % Initial conditions
dvds = @(s,v_vec)state_deriv(s,v_vec,param);
ds = 1/(N-1);
[s, v_mat] = rk4(dvds, v0, ds, N-1); % Time step state system

% Grid-stretching ratio function
f1 = (A/2)*(1-erf(b*(s-c)));
fGSR = f1 - A;

% Physical domain
y = v_mat(:,1);
dy = y(2:length(y)) - y(1:length(y)-1);
gsr = (dy(2:length(dy)) - dy(1:length(dy)-1))./dy(1:length(dy)-1);
fprintf("fGSR mesh results for N = %d:\n",length(y));
fprintf("dy/ds(s=0)= %6.3E\n",vfinal);
fprintf("dy max    = %6.3E\n",max(dy));
fprintf("dy min    = %6.3E\n",min(dy));
fprintf("max/min   = %6.3E\n",max(dy)/min(dy));
fprintf("max gsr   = %6.3f%%\n",max(gsr)*100);

% plot results
figure() % plot fGSR(s)
subplot(2,1,1) 
plot(s,fGSR,'-k');
str = sprintf("fGSR Y mesh; N = %d",N);
xlabel('Computational domain [s]');
ylabel('Grid-stretching ratio [fGSR]');
title(str);
subplot(2,1,2) 
plot(s,y,'or');
xlabel('Computational domain [s]');
ylabel('Physical domain [y]');
title(str);
visualize_xy_mesh(0,1,ymin,ymax,[0,1],y);

%% Problem 3
kdx = linspace(0,pi,200);           % common domain
exact = kdx;                        % exact solution
kdx3re = sin(kdx).*(2-cos(kdx));    % real component
kdx3im = sin(kdx).^2+2*cos(kdx)-2;  % imaginary

figure()
plot(kdx,exact,'-k');
hold on
plot(kdx,kdx3re,'--r');
plot(kdx,kdx3im,'--b');
xlim([0 pi]);
ylim([-pi pi]);
legend('Exact','Upwind Real','Upwind Im');
title('Upwind approximation');
xlabel('$$k\Delta x$$','Interpreter','latex');
ylabel('$$\tilde{k}\Delta x$$','Interpreter','latex');

%% Functions
% Derivative function for state system
% Parameters:
%   s     - location value to evaluate dvds at
%   v_vec - state vector to evaluate dvds at
% Outputs:
%   dvds  - the output vector of derivatives of state vector
% Calling:
%   @(s,v_vec)state_deriv(s, v_vec); params are passed in
function dvds = state_deriv(s, v_vec,param)
    % parameters passed in
    Ny = param(3);  % Number of y points
    A = param(4);   % Tuning parameters
    b = param(5);   %
    c = param(6);   %
    
    f1 = (A/2)*(1-erf(b*(s-c)));
    fGSR = f1 - A;             % Grid-stretching function
     
    n = size(v_vec,1); % v_vec comes in as column state vector
    dvds = zeros(n,1); % dvds will leave as column state vector
    dvds(1) = v_vec(2);
    dvds(2) = (Ny-1)*fGSR*v_vec(2);
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

% Computes the residual of grid-stretching second-order DE for vguess
% Parameters:
%   vguess - independent variable for shooting parameter
% Calling:
%   residual expected to be called as R(vguess); params are passed in
function res = residual(vguess,param)
    % parameters passed in
    ymin = param(1);        % initial condition [m]
    N = param(3);           % Number of points
    ymax = param(2);        % Boundary Condition for y at s=1
    
    % Run rk4 to time-step state vector to s=1 using vguess
    v0 = [ymin;vguess];    % Initial conditions
    dvds = @(s,v_vec)state_deriv(s,v_vec,param);
    ds = 1/(N-1);
    [s, v_mat] = rk4(dvds, v0, ds, N-1); % Time step state system
    
    % Pull out end results from rk4 (at edge of computational domain s=1)
    yend = v_mat(end,1);
   
    % Compute residual
    res = yend - ymax;
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

% Runge-Kutta approximation with 4 slopes
% Parameters:
%   dv_fcn - function handle which evaluates derivative of state vector
%   v0_fcn - vector with initial value of each state variable [nx1]
%   dt     - time step 
%   N      - number of time steps
% Outputs:
%   t_vec  - vector of time values [ (N+1) x 1 ]
%   v_mat  - matrix with instantaneous value of each state variable 
%            [ (N+1) x n ]
% Calling:
%   dv_fcn expected to be called as dv_fcn(t, v_vec), where v_vec is the
%   current state vector
function [t_vec, v_mat] = rk4(dv_fcn, v0_vec, dt, N)
    n = size(v0_vec,1);     % n = number of state variables
    t_vec = (0:dt:dt*N)';   % time column vector
    v_mat = zeros(N+1,n);   % columns for each state var; rows for time steps
    v_mat(1,:) = v0_vec';   % row=0; t=0; v = v0 (ICs)
    for i = 1:N
        k1 = dv_fcn(t_vec(i),v_mat(i,:)'); % d{v}/dt 1 (col vector)
        ystar1 = v_mat(i,:)' + (dt/2)*k1;
        
        k2 = dv_fcn(t_vec(i) + dt/2, ystar1); % d{v}/dt 2
        ystar2 = v_mat(i,:)' + (dt/2)*k2;
        
        k3 = dv_fcn(t_vec(i) + dt/2, ystar2); % d{v}/dt 3
        ystar3 = v_mat(i,:)' + dt*k3;
        
        k4 = dv_fcn(t_vec(i) + dt, ystar3); % d{v}/dt 4
        v_mat(i+1,:) = v_mat(i,:) + dt*(k1' + 2*k2' + 2*k3' + k4')/6;
    end
end

% Plot an xy-plane structured mesh
% Code provided by Dr. Lui
% Parameters:
%   x_min - minimum x boundary
%   x_max - maximum y boundary
%   y_min - minimum y boundary
%   y_max - maximum y boundary
%   x_vec - vector containing all x mesh points
%   y_vec - vector containing all y mesh points
% Outputs:
%   Plots the resulting mesh on the xy-plane
function [] = visualize_xy_mesh(x_min, x_max, y_min, y_max, x_vec, y_vec)

figure()
hold on;
Nx = length(x_vec);
Ny = length(y_vec);

for j = 1:Ny
    x_horizontal(j,:) = [x_min x_max];
    y_horizontal(j,:) = [y_vec(j) y_vec(j)];
    plot(x_horizontal(j,:), y_horizontal(j,:));
end

for i = 1:Nx
    x_vertical(i,:) = [x_vec(i) x_vec(i)];
    y_vertical(i,:) = [y_min y_max];
    plot(x_vertical(i,:), y_vertical(i,:));
end

hold off;
axis([x_min x_max y_max y_min]);
xlabel('x');
ylabel('y');
title(sprintf('Structured Mesh for a Free Jet (Ny = %3i)',Ny));
end