% initialization
clear all
close all
clc
% data section
tf = 3;        % s
mp = .0685;    % kg
mw = .088;    % kg
Lp = .432;   % m
dw = .05;   % m
ds = .025;  % m
g = 9.8;    % m/s^2
Lpcg = (Lp - ds)/2;   % m
tol = 1*10^(-6);     % See Simulink->Simulation-> 
maxstep = 0.01;      % Model Configuration Parameters->Additional Options

% a(t) case 1
T=[0:0.01:tf]';   % s
x = -43.5668.*T + 9.38;
a1 = (7194.*exp(x).*(exp(x)-1))./(1+exp(x)).^3;  % cm/s^2

% a(t) case 2
y = -43.5668*T + 31.1633;
a2 = (3597.*exp(x).*(exp(x)-1))./(1+exp(x)).^3 + (3597.*exp(y).*(exp(y)-1))./(1+exp(y)).^3;

% vectors
maxVelocitiesA1 = 0.04:0.001:0.39;   % using a(t) case 1
maxVelocitiesA2 = 0.04:0.001:0.39;   % using a(t) case 2
offsets = 0.04:0.001:.39;            % Lwcg offset
n = 1;

% finds max angular velocity of pendulum's residual swing time as 
% a function of weight offset distance (cm) and stores results in vectors
% (for each a(t) case)
for Lwcg=.04:0.001:.39
    J = mp*Lp^2/12 + mp*Lpcg^2 + (1/2)*mw*(dw/2)^2 + mw*Lwcg^2;
    k = g*(mp*Lpcg+mw*Lwcg);
    C = mp*Lpcg + mw*Lwcg;

    sim('crane_simulink_a1')
    offsets(n) = 100 * Lwcg;   % cm
    maxVelocitiesA1(n) = max(omega(find(t>1)));
    
    sim('crane_simulink_a2')
    maxVelocitiesA2(n) = max(omega(find(t>1)));
    
    n = n + 1;
end
% find minimum value
[minVelocity, index] = min(maxVelocitiesA2);
minoffset = offsets(index);

% verify result
fprintf('Min velocity is %6.4f. Corresponding offset is %6.4f.', minVelocity, minoffset);
Lwcg = minoffset;
sim('crane_simulink_a2')
figure
plot(t, omega);
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')

save Lab7MaxResidualAngularVelocityA1.mat offsets maxVelocitiesA1;
save Lab7MaxResidualAngularVelocityA2.mat offsets maxVelocitiesA2;

% offset v. max residual angular velocity (a1)
figure
plot(offsets, maxVelocitiesA1);
xlabel('Moveable weight offset (cm)');
ylabel('Residual swing maximum angular velocity (rad/s)');

% offset v. max residual angular velocity (a2)
figure
plot(offsets, maxVelocitiesA2);
xlabel('Moveable weight offset (cm)');
ylabel('Residual swing maximum angular velocity (rad/s)');

% other (optional) plots

% sim('crane_simulink')
% figure
% plot(t, omega);
% xlabel('Time (s)')
% ylabel('Angular Velocity (rad/s)')
