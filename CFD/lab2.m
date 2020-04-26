%% lab2.m
% Joshua Eckels
% ME427 - Intro to CFD
% Lab 2 - Meshing with grid-spacing fGS and grid-stretching ratio fGSR
% 3/30/20
clc
close all
clear all
% Perform the grid-meshing for these N values
N = [25,50,100,200];
fig_idx = 1;
for i = 1:length(N)
    %% Streamwise mesh (x-direction) using fGS for grid-spacing
    % Tuning parameters for fGS = A*exp(alpha*s)
    alpha = 2;
    xmin = 0;
    xmax = 40;
    A = xmax*alpha/(exp(alpha)-1);
    
    % computational domain
    ds = 1/(N(i)-1);            % domain spacing
    s = 0:ds:1;                 % domain range
    fGS = A*exp(alpha*s);       % grid-spacing function fGS
    
    % physical domain
    x = (A/alpha)*exp(alpha*s)-(A/alpha);
    dx = x(2:length(x)) - x(1:length(x)-1);
    gsr = (dx(2:length(dx)) - dx(1:length(dx)-1))./dx(1:length(dx)-1);
    fprintf("Streamwise (x) mesh results for N = %d:\n",length(x));
    fprintf("dx max   = %4.3f\n",max(dx));
    fprintf("dx min   = %4.3f\n",min(dx));
    fprintf("max/min  = %4.3f\n",max(dx)/min(dx));
    fprintf("max gsr  = %4.3f%%\n",max(gsr)*100);
    
    % Analytical dx/ds and d2x/ds2
    dxds = A*exp(alpha*s);
    d2xds2 = A*alpha*exp(alpha*s);
    
    % plot results
    figure(1) % plot fGS(s)
    subplot(2,2,i) 
    plot(s,fGS,'.k');
    str = sprintf("Streamwise mesh; N = %d",N(i));
    xlabel('Computational domain [s]');
    ylabel('Grid-spacing [fGS]');
    title(str);
    figure(2) % plot x(s)
    subplot(2,2,i) 
    plot(s,x,'.k');
    xlabel('Computational domain [s]');
    ylabel('Physical domain [x]');
    title(str);
    figure(3) % plot dx/ds(s)
    subplot(2,2,i) 
    plot(s,dxds,'.k');
    xlabel('Computational domain [s]');
    ylabel('Derivative [dx/ds]');
    title(str);
    figure(4) % plot d^2x/ds^2(s)
    subplot(2,2,i) 
    plot(s,d2xds2,'.k');
    xlabel('Computational domain [s]');
    ylabel('Second derivative [d2x/ds2]');
    title(str);
    
    %% Transverse mesh (y-direction) using fGSR for grid-stretching ratio
    % parameters
    ymin = 0;                               % Lower y bound [m]
    A = 1/4;                                % fGSR tuning parameters
    b = 4;                                 %
    c1 = 1/8;                               %
    c2 = 7/8;                               %
    ymax = 1;                               % Upper y bound [m]
    param = [ymin,N(i),A,b,c1,c2,ymax];

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
    ds = 1/(N(i)-1);
    [s, v_mat] = rk4(dvds, v0, ds, N(i)-1); % Time step state system
    
    % Grid-stretching ratio function
    f1 = (A/2)*(1-erf(b*(s-c1)));
    f2 = (A/2)*(1-erf(b*(s-c2)));
    fGSR = f1 + f2 - A;
    
    % Physical domain
    y = v_mat(:,1);
    dyds = v_mat(:,2);
    d2yds2 = (N(i)-1)*fGSR.*dyds;
    dy = y(2:length(y)) - y(1:length(y)-1);
    gsr = (dy(2:length(dy)) - dy(1:length(dy)-1))./dy(1:length(dy)-1);
    fprintf("Transverse (y) mesh results for N = %d:\n",length(y));
    fprintf("dy/ds(s=0)= %6.3E\n",vfinal);
    fprintf("dy max    = %6.3E\n",max(dy));
    fprintf("dy min    = %6.3E\n",min(dy));
    fprintf("max/min   = %6.3E\n",max(dy)/min(dy));
    fprintf("max gsr   = %6.3f%%\n",max(gsr)*100);
    
    % plot results
    figure(5) % plot fGSR(s)
    subplot(2,2,i) 
    plot(s,fGSR,'.k');
    str = sprintf("Transverse mesh; N = %d",N(i));
    xlabel('Computational domain [s]');
    ylabel('Grid-stretching ratio [fGSR]');
    title(str);
    figure(6) % plot y(s)
    subplot(2,2,i) 
    plot(s,y,'.k');
    xlabel('Computational domain [s]');
    ylabel('Physical domain [y]');
    title(str);
    figure(7) % plot dy/ds(s)
    subplot(2,2,i) 
    plot(s,dyds,'.k');
    xlabel('Computational domain [s]');
    ylabel('Derivative [dy/ds]');
    title(str);
    figure(8) % plot d^2y/ds^2(s)
    subplot(2,2,i) 
    plot(s,d2yds2,'.k');
    xlabel('Computational domain [s]');
    ylabel('Second derivative [d2y/ds2]');
    title(str);
    
    % Visualize mesh in xy-plane (for N=25)
    if (i == 1)
        visualize_xy_mesh(xmin,xmax,ymin,ymax,x,y);
    end
end


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
    Ny = param(2);  % Number of y points
    A = param(3);   % Tuning parameters
    b = param(4);   %
    c1 = param(5);  %
    c2 = param(6);  %
    
    f1 = (A/2)*(1-erf(b*(s-c1)));
    f2 = (A/2)*(1-erf(b*(s-c2)));
    fGSR = f1 + f2 - A;             % Grid-stretching function
     
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
    N = param(2);           % Number of points
    ymax = param(7);        % Boundary Condition for y at s=1
    
    % Run rk4 to time-step state vector to s=1 using vguess
    v0 = [ymin;vguess];    % Initial conditions
    dvds = @(s,v_vec)state_deriv(s,v_vec,param);
    ds = 1/(N-1);
    [s, v_mat] = rk4(dvds, v0, ds, N-1); % Time step state system
    
    % Pull out end results from rk4 (at edge of computational domain s=1)
    yend = v_mat(end,1);
    dydsend = v_mat(end,2);
   
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

figure(9) % figure number hard-coded
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
axis([x_min x_max y_min y_max]);
xlabel('x');
ylabel('y');
title(sprintf('Structured Mesh for a Free Jet (Nx = %3i, Ny = %3i)', Nx, Ny));
end