%% Eckels_Joshua_L9.m
% ME327
% Joshua Eckels
% 2/14/2019
clear all
close all
clc
%% Task 2
v0 = 2;     % IC
dt = 0.1;   % time step
dvdt = @(t,v_vec) 2*(t - v_vec);
[t_vec, v_mat] = rk4(dvdt,v0,dt,10);
fprintf("Task 2\n");
fprintf("The Runge-Kutta method result for Task 2 IVP at t=0.1 is v=%6.4f\n",v_mat(2));

%% Task 3
dvdt = @(t,v_vec) 2*(sin(t) - v_vec);
v0 = 2;
dt = [0.8;0.08;0.008];
tf = 8;
v_exact = 2.4*exp(-2*tf)+0.8*sin(tf)-0.4*cos(tf);
v_approx = zeros(3,1);
error = zeros(3,1);
fprintf("\nTask 3:\n");
fprintf("%10s %10s %14s\n","Time step","y(t = 8)","error"); 
for i = 1:3
    N = tf/dt(i);  % number of time steps
    [t_vec, v_mat] = rk4(dvdt,v0,dt(i),N);
    v_approx(i) = v_mat(end,1);
    error(i) = v_exact - v_approx(i);
    fprintf("%10.4f %10.8f %14.11f\n",dt(i),v_approx(i),error(i));
end
figure(1)
plot(t_vec,v_mat(:,1),'.k');
hold on
t = linspace(0,tf,1000);
v = 2.4*exp(-2*t)+0.8*sin(t)-0.4*cos(t);
plot(t,v,'-r');
legend("RK4","Analytic");
title("RK4 approximation");
xlabel("Time t [s]");
ylabel("Displacement v");

%% Task 5
v0 = [pi/2;0];      % Initial conditions for pendulum
dvdt = @(t,v_vec)pend_deriv(t,v_vec);
tf = 3;             % simulate for 3 seconds
dt = 0.01;          % time step of 0.01 seconds
N = tf/dt;          % number of time steps
[t_vec, v_mat] = rk4(dvdt,v0,dt,N);

figure(2)
title('RK4 Pendulum approximation');
hold on
grid on
theta = v_mat(:,1);
omega = v_mat(:,2);
subplot(211);
plot(t_vec,theta,'.k');
xlabel('Time [s]');
ylabel('Theta [rad]');
subplot(212);
plot(t_vec,omega,'.k');
xlabel('Time [s]');
ylabel('Theta dot [rad/s]');

fprintf("\nTask 5\n");
fprintf("At t=%d, theta = %8.6f rad, omega = %8.6f rad/s\n",t_vec(end),theta(end),omega(end));

%% Task 6
% approximate the error with Richardson's Extrapolation
% result from N=300 saved above, repeat with N=600
dt2 = dt/2;         % divide time step in half
N2 = tf/dt2;        % number of time steps
[t_vec2, v_mat2] = rk4(dvdt,v0,dt2,N2);

% find error from N=300 to N=600 for both theta and omega
theta_diff = zeros(1,size(t_vec,1));
omega_diff = zeros(1,size(t_vec,1));
n = dt/dt2; % amount dt was divided by from N=300 to N=600
for i = 1:size(v_mat,1)
    theta_diff(i) = v_mat2(((i-1)*2+1),1) - v_mat(i,1);
    omega_diff(i) = v_mat2(((i-1)*2+1),2) - v_mat(i,2);
end
theta_error = (n^4/(n^4-1))*(theta_diff);
omega_error = (n^4/(n^4-1))*(omega_diff);
theta_norminf = max(abs(theta_error));
omega_norminf = max(abs(omega_error));
norminf_avg = (theta_norminf + omega_norminf) / 2;

fprintf("\nTask 6\n");
fprintf("Error in pendulum angle at t = %d s is estimated to be %6.4e rad.\n", t_vec(end),theta_norminf);
fprintf("Error in pendulum velocity at t = %d s is estimated to be %6.4e rad/s.\n", t_vec(end),omega_norminf);
fprintf("The average error at t = %d s is estimated to be %6.4e.\n", t_vec(end),norminf_avg);

%% Task 7
options = odeset('AbsTol',theta_norminf);
[t_ode, v_ode] = ode45(dvdt, [0 tf], v0,options);

figure(3)
title('RK4 v. ode45 Comparison');
grid on
subplot(211);
plot(t_vec,theta,'-r');
hold on
plot(t_ode,v_ode(:,1),'.k');
xlabel('Time [s]');
ylabel('Theta [rad]');
legend('RK4', 'ode45');
subplot(212);
plot(t_vec,omega,'-r');
hold on
plot(t_ode,v_ode(:,2),'.k');
xlabel('Time [s]');
ylabel('Theta dot [rad/s]');
legend('RK4', 'ode45');

fprintf("\nTask 7\n");
fprintf("RK4 number of time steps is %d. Ode45 number of time steps is %d\n",size(t_vec,1),size(t_ode,1));

figure(4)
temp = t_ode(2:end,1);
tstep = temp - t_ode(1:end-1,1);
indices = (1:size(tstep,1))';
semilogy(indices,tstep,'ok');
xlabel('Time-step index');
ylabel('Time step (dt)')

%% Functions

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

% Derivative function for pendulum problem in Task 4
% Parameters:
%   t     - time value to evaluate dvdt at
%   v_vec - state vector to evaluate dvdt at
% Outputs:
%   dvdt  - the output vector of derivatives of state vector
% Calling:
%   @(t,v_vec)pend_deriv(t, v_vec)
function dvdt = pend_deriv(t, v_vec)
     mp = 0.0685;       % mass of rod [kg]
     Lp = 0.432;        % length of rod [m]
     mw = 0.088;        % mass of weight [kg]
     dw = 0.05;         % dia of weight [m]
     ds = 0.025;        % dia of sensor [m]
     Lwcg = 0.30;       % location of weight [m]
     
     Lpcg = (Lp-ds)/2;  % center of mass of pole [m]
     k = 9.8*(mp*Lpcg+mw*Lwcg);
     J = mp*Lp^2/12+mp*Lpcg^2+0.5*mw*(dw/2)^2+mw*Lwcg^2;
     
     n = size(v_vec,1); % v_vec comes in as column state vector
     dvdt = zeros(n,1); % dvdt will leave as column state vector
     dvdt(1) = v_vec(2);
     dvdt(2) = (-k/J)*sin(v_vec(1));
end