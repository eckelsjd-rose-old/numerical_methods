%% lab5.m
% ME327 - Numerical Methods
% Joshua Eckels
clear all
close all
clc

% search interval
r_lower = 0;
r_upper = 0.01
h_lower = 0;
h_upper = 0.04;

% known values
wn_allow = 100;          % natural frequency [rad/s]
% orignal value of 50 was too small
ts_allow = 2;            % settling time [s]
rho = 8000;              % density [kg/m^3]
Cm = 0.1*1000;           % cost per kg [$/kg]
k = 500;                 % spring constant [N/m]
c = 0.5;                 % damping coefficient [N-s/m]

% design variables
rVec = (r_lower:0.001:r_upper)';    % radius [m]
hVec = (h_lower:0.001:h_upper)';    % height [m]

% set optimal value to an unreasonably high value
f_min = 1e15;

% initialize matrices
rMat = zeros(size(rVec,1),size(hVec,1));
hMat = zeros(size(rVec,1),size(hVec,1));
fMat = zeros(size(rVec,1),size(hVec,1));

tic
% compute the objective function
for i = 1:length(rVec) 
    r = rVec(i);
    for j = 1:length(hVec)
        h = hVec(j);
        
        rMat(i,j) = r;
        hMat(i,j) = h;
        
        % check constraints
        wn = sqrt(k/(rho*pi*r^2*h));
        ts = 8/(c/(rho*pi*r^2*h));
        if ((ts-ts_allow <= 0) && (wn-wn_allow <= 0))
            fMat(i,j) = Cm*rho*pi*r^2*h;
        else
            fMat(i,j) = NaN;
        end
        
        if fMat(i,j) < f_min
            f_min = fMat(i,j);
            r_min = r;
            h_min = h;
        end
    end
end
toc
    
% plot result as a surface
figure(1)
clf
surf(rMat, hMat, fMat, 'EdgeAlpha',0);
xlabel('$r$','Interpreter','Latex');
ylabel('$h$','Interpreter','Latex');
zlabel('$f(r,h)$','Interpreter','Latex');

fprintf("The minimum cost of $%6.4f is achieved when r=%6.3f m and h=%6.3f m.\n",f_min,r_min,h_min);