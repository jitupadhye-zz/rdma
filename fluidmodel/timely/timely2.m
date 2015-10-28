clc;clear all;close all;

global C; % bandwidth 
global Seg; % MSS 
global delta; % the additive increment step
global T_high; % if RTT is greater than this, decrease rate multiplicatively.
global T_low; % if RTT is lower than this, do increase rate additively. 
global prop; % propagation delay. 
global minRTT; % 20 microseconds, defined protocol parameter
global beta; % beta, protocol parameter
global alpha; % alpha, protocol parameter.
global callCounter;
global maxQueue;

%
% Just for testing purposes.
%
callCounter = 0;

% 
% Fixed Parameters
%
step_len = 5/1e6 ; % 5 microseconds
sim_length = step_len * 10000 ;
C = 10 * 1e9; % 10Gbps
Seg = 16 * 8 * 1e3; % 16KB
prop = 5 /1e6; % 5 microseconds - picked double the lowest one way delay we have seen. 

%
% Parameters we can play with.
%
delta = 10 * 1e6; % 10Mbps, expressed in Kilobytes per second.
T_high = 500/1e6; % 500 microseconds (see section 4.4)
T_low = 50/1e6; % 50 microseconds (see section 4.4). 
minRTT = 20/1e6; % 20 microseconds 
beta = 0.8;
alpha = 0.875; % unsure
maxQueue = 2 * C * T_high; % only for corner cases - queue won't grow beyond this. 

%
% Initial conditions: 
%
ir1 = 19 * 1e9; % initial rate of flow 1
ir2 = 0 * 1e9; %  initial rate of flow 2
q0 =  0; % initial queue length

%
% Options.
%
options = ddeset('MaxStep', step_len);

%
% Solve.
%
sol = ddesd(@fluid_timely, @fluid_delays, [ir1; 0; ir2; 0; q0], [0, sim_length],options);

%
% Extract solution.
%
t = sol.x;
r = sol.y(1,:);
g = sol.y(2,:);
r2 = sol.y(3,:);
g2 = sol.y(4,:);
q = sol.y(5,:);

%
% Write solution to file.
%
dlmwrite('timely_fluid.txt',[t',r'./1e9, r2'./1e9, q'./(8*1e3)],'\t');

%
% Plot.
%
figure

subplot(1,2,1);
plot(t,r./1e9,'b');
hold on
plot(t,r2./1e9,'r--');
%axis([0 sim_length 0 max(max(r./1e9), max(r2./1e9))])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
xlabel('Time (seconds)')
ylabel('Throughput (Gbps)')

%plot queue size
subplot(1,2,2);
plot(t,q./(8*1e3))
xlabel('Time (seconds)')
ylabel('Queue (KBytes)')

set(gcf,'Units','centimeters','Position',[0 0 35 14])
