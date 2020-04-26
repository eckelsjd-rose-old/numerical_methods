%% lab4.m
% Class:       ME427 - Intro to CFD; Lab 4
% Author:      Joshua Eckels
% Description: Function to evaluate non-linear, 1D, steady Burgers
%              equation with a non-uniform mesh
% Strategy:    1. Compute a non-uniform mesh
%              2. Define the finite volume conservation equation
%              3. Use interpolation or Taylor to get interface values
%              4. Define system of N residual eqns for N unknown grid cells
%              5. Solve system with multi-dimensional Newton-Raphson method
% Date:        4/9/2020
clc
clear all
close all
%% Solve the Burgers equation for u(x), with nu = 0.1
N1 = 150;
nu = 0.1;
A = -0.05;
b = 30;
c = 0.5;
fprintf("TASK 1:\n");
fprintf("Finite volume solution with nu=%3.1f\n",nu);

% generate the mesh
[s,x_mat,fGSR] = SOLVE.generate_mesh(N1+1,A,b,c);
figure() % plot fGSR(s)
subplot(2,1,1) 
plot(s,fGSR,'-r');
xlabel('Computational domain [s]');
ylabel('Grid-stretching ratio [fGSR]');
title(sprintf("X mesh for N = %d",N1));
subplot(2,1,2) 
plot(s,x_mat(:,1),'-r');
xlabel('Computational domain [s]');
ylabel('Physical domain [x]');
SOLVE.display_mesh(x_mat(1,1), x_mat(end,1), -1, 1, x_mat(:,1), [-1, 1]);

[x1, urt1] = SOLVE.solve_burgers(N1,nu,A,b,c);
dx = x1(2:length(x1)) - x1(1:length(x1)-1);
x_cells = x1(1:length(x1)-1) + dx/2;

% load finite-difference result for nu = 0.1
fd_data = load('fd.mat');
x_fd = fd_data.x1;
u_fd = fd_data.urt1;

figure()
plot(x_cells,urt1,'or');    % finite volume solution
hold on
plot(x_fd,u_fd,'-k');       % finite difference solution
xlabel('X coordinate [m]');
ylabel("Function u [m/s]");
title(sprintf('Solution to 1D Burgers equation; N=%d, nu=%3.1f',N1,nu));
xlim([-1 1]);
legend('FV','FD');
grid on

% interpolate FV and FD solutions to common domain
N2 = size(x_fd,1);
x_common = linspace(x1(1),x1(end),(N1+N2)/2);
u1c = interp1(x_cells,urt1,x_common);
u2c = interp1(x_fd,u_fd,x_common);

% approximate error in N1 finite volume solution
diff = u1c - u2c;
norminf = max(abs(diff));
fprintf("Error at N=%d is estimated to be %6.3E m/s.\n",N1,norminf);

%% Solve the Burgers equation for u(x), with nu = 0.0001
N1 = 180;
nu = 0.0001;
A = -0.09;
b = 30;
c = 0.5;
fprintf("\nTASK 2:\n");
fprintf("Finite volume solution with nu=%3.4f\n",nu);

% generate the mesh
[s,x_mat,fGSR] = SOLVE.generate_mesh(N1+1,A,b,c);
figure() % plot fGSR(s)
subplot(2,1,1) 
plot(s,fGSR,'-r');
xlabel('Computational domain [s]');
ylabel('Grid-stretching ratio [fGSR]');
title(sprintf("X mesh for N = %d",N1));
subplot(2,1,2) 
plot(s,x_mat(:,1),'-r');
xlabel('Computational domain [s]');
ylabel('Physical domain [x]');
SOLVE.display_mesh(x_mat(1,1), x_mat(end,1), -1, 1, x_mat(:,1), [-1, 1]);

[x1, urt1] = SOLVE.solve_burgers(N1,nu,A,b,c);
dx = x1(2:length(x1)) - x1(1:length(x1)-1);
x_cells = x1(1:length(x1)-1) + dx/2;

% load finite-difference result for nu = 0.0001
fd_data = load('fd2.mat');
x_fd = fd_data.x1;
u_fd = fd_data.urt1;

figure()
plot(x_cells,urt1,'or');    % finite volume solution
hold on
plot(x_fd,u_fd,'-k');       % finite difference solution
xlabel('X coordinate [m]');
ylabel("Function u [m/s]");
title(sprintf('Solution to 1D Burgers equation; N=%d, nu=%3.4f',N1,nu));
xlim([-1 1]);
legend('FV','FD');
grid on

% interpolate FV and FD solutions to common domain
N2 = size(x_fd,1);
x_common = linspace(x1(1),x1(end),(N1+N2)/2);
u1c = interp1(x_cells,urt1,x_common);
u2c = interp1(x_fd,u_fd,x_common);

% approximate error in N1 finite volume solution
diff = u1c - u2c;
norminf = max(abs(diff));
fprintf("Error at N=%d is estimated to be %6.3E m/s.\n",N1,norminf);