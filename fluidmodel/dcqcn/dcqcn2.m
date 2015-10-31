clc;clear all;close all;

global Rai;
global C;
global B;
global F;
global g;
global Kmin;
global Kmax;
global timer;
global pmax;
global tau;
global tauprime;

%
% simulation control.
%
sim_step = 5 * 1e-6; % 5 microseconds.
sim_length = sim_step * 20000; 

%
% Fixed parameters. 
%
C = 40 * 1e9;   % 40Gbps. Link speed. 
tau = 50 * 1e-6;  % 50 microseconds. This is the feedback delay. 
tauprime = 55 * 1e-6; % 55 microseconds. This is the interval of equation 2. 
F = 5; % Fast recovery steps. 
B = 10 * 8 * 1e6;   %10MB.Byte counter.
Rai = 40 * 1e6; % 40Mbps. Rate increase step.

%
% Tunable parameters.
%
timer = 55 * 1e-6; % 55 microseconds. This is rate increase timer.
Kmax = 200 * 8 * 1e3; % 200KB
Kmin = 5 * 8 * 1e3; % 5KB
pmax = 1e-2; % 1 percent. 
g = 1/256; 


%
% Initial conditions. 
%
ir1 = C;
ir2 = 0;

%
% solve 
%
options = ddeset('MaxStep', sim_step);
sol = dde23('fluid2_test', tau, [ir1; C; 1; ir2; C; 1; 0], [0, sim_length], options);

%
% parse output.
%
t = sol.x;
rc = sol.y(1,:);
rt = sol.y(2,:);
alpha = sol.y(3,:);
rc2 = sol.y(4,:);
rt2 = sol.y(5,:);
alpha2 = sol.y(6,:);
q = sol.y(7,:);

%
% write to file.
%
dlmwrite(sprintf('fluid.txt'),[t',rc'./1e9,rc2'./1e9, q'./(8*1e3)],'\t');

% 
% plot.
%
figure
subplot(1,2,1);
plot(t, rc./1e9, 'b', t, rc2./1e9, 'r--');
xlabel('Time (seconds)')
ylabel('Throughput (Gbps)')
subplot(1,2,2);
plot(t,q./(8*1e3))
xlabel('Time (seconds)')
ylabel('Queue (KBytes)')

