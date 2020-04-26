%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% me411_dp2_script.m
%
%  THIS PROGRAM
%  This program iterates through several combinations of bypass ratios
%  and compressor ratios and generates a carpet plot of specific thrust
%  v. thrust-specific fuel consumption for a turbofan engine.
%
% Input: All input is in the initialization section.
% Output: Prints to a carpet plot.
%
%   Written by Joshua Eckels
%              4/4/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % initialization % % %
clear all
close all
clc

% % % % % % % % % % % % % % % given values % % %
T0 = 216.775;               % K
M0 = 0.83;                  %
hPR = 42800;                % kJ/kg
Cp = 1.004;                 % kJ/kg-K
R = 287;                    % J/kg-K
g = 1.4;                    %
Tt4 = 1560;                 % K
prF = 1.5;                  % fan pressure ratio
prC = 20:4:40;              % compressor pressure ratio

astar = 1;                  % initialize optimum bypass ratio
tauC = 1;                   % initialize compressor temp ratio
alpha=[4,6,8,10,12,astar];  % bypass ratios

% % % % % % % % % % % % % % % changing values % % %
Tt3 = zeros(6, 6);          % Compressor outlet temperature
f = zeros(6, 6);            % mf/mc
Tt5 = zeros(6, 6);          % turbine outlet temperature
prT = zeros(6, 6);          % turbine pressure ratio
pr9 = zeros(6, 6);          % core exit pressure ratio
M9 = zeros(6, 6);           % core exit Mach
T9 = zeros(6, 6);           % core exit temp
v9 = zeros(6, 6);           % core exit velocity
Fm0 = zeros(6, 6);          % Specific Thrust F/m0
S = zeros(6, 6);            % Thrust specific fuel consumption

% % % % % % % % % % % % % % % constant values % % %
v0 = M0*sqrt(g*R*T0);       % Inlet velocity
prR = (1 + 0.2*M0^2)^3.5;   % Inlet pressure ratio
Tt0 = (1 + 0.2*M0^2)*T0;    % Inlet temperature
Tt13 = prF^(1/3.5)*Tt0;     % Exit fan temperature
pr19 = prF*prR;             % Fan nozzle pressure ratio
M19 = sqrt((pr19^(1/3.5)-1)/0.2);
T19 = Tt13/(1+0.2*M19^2);   % Fan exit temp
v19 = M19*sqrt(g*R*T19);    % Fan exit velocity
tauR = 1 + 0.2*M0^2;        % Inlet temp ratio
tauL = Tt4/T0;              % Turbine temp ratio
tauF = prF^(1/3.5);         % Fan temp ratio

% % % % % % % % % % % % % % % iterate alpha and prC % % %
for i = 1:length(prC)
    for j = 1:length(alpha)
        % calculate optimum bypass ratio
        tauC = prC(i)^(1/3.5);
        alpha(6) = (1/(tauR*(tauF-1)))*(tauL-tauR*(tauC-1)-tauL/(tauR*tauC)-(1/4)*(sqrt(tauR*tauF-1)+sqrt(tauR-1))^2);
        % calculate intermediate values
        Tt3(i,j) = prC(i)^(1/3.5)*Tt0;
        f(i,j) = Cp*(Tt4-Tt3(i,j))/hPR;
        Tt5(i,j) = -1*(alpha(j)*(Tt13-Tt0)+Tt3(i,j)-Tt0-Tt4);
        prT(i,j) = (Tt5(i,j)/Tt4)^3.5;
        pr9(i,j) = prT(i,j)*prC(i)*prR;
        M9(i,j) = sqrt((pr9(i,j)^(1/3.5)-1)/0.2);
        T9(i,j) = Tt5(i,j)/(1+0.2*M9(i,j)^2);
        v9(i,j) = M9(i,j)*sqrt(R*g*T9(i,j));
        % outputs
        Fm0(i,j) = (v9(i,j)-v0 + alpha(j)*(v19-v0))/(1+alpha(j)); % (N/(kg/s))
        S(i,j) = f(i,j)*10^6/((1+alpha(j))*Fm0(i,j)); % ((g/s)/kN)
    end    
end

% % % % % % % % % % % % % % % carpet plot % % %
figure;
hold on;
box on;
for m = 1:length(prC)
    for n = 1:length(alpha)
        plot(Fm0(m,n),S(m,n),'or');
        plot(Fm0(:,n),S(:,n),'r-');
    end
    plot(Fm0(m,:),S(m,:),'k-');
end
% plot aircraft performance limitations
Smax = 16.5424;             % (g/s)/kN
Fm0min = 139.63;            % N/(kg/s)
line([Fm0min, Fm0min], get(gca, 'ylim'));
line(get(gca, 'xlim'), [Smax,Smax]);

xlabel('F/m0 [N/(kg/s)]');
ylabel('S [(g/s)/kN]');
