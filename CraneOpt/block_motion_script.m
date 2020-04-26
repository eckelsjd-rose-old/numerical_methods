clear all
close all
clc
m = 2;    % kg
tf = 4;   % s
v0 = 1;   % m/s

T=[0:0.01:tf]';
f=10.*exp(-0.6.*T).*sin(8.*T);   % N

sim('block_motion_sim')

figure
plot(t,v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')