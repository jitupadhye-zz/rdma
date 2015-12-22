function timely2()

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
global maxQueue; % max queue. 
global numFlows;  % number of flows. 
global initVal; % column of initial values (rate and RTT gradient for each flow, plus initial queue length). 

%
% Simulation control
% 
step_len = 5e-6 ; % 5 microseconds
sim_length = step_len * 20000 ;


% 
% Fixed Parameters
%
C = 10 * 1e9; % line rate.
Seg = 64 * 8 * 1e3; % burstsize.
prop = 4e-6; % propagation delay
numFlows = 3; % number of flows. 

%
% Parameters we can play with.
%
delta = 10e6; % 10Mbps
T_high = 500e-6; % 500 microseconds (see section 4.4)
T_low = 50e-6; % 50 microseconds (see section 4.4). 
minRTT = 20e-6; % 20 microseconds 
beta = 0.8;
alpha = 0.875; % unsure
maxQueue = 2 * C * T_high; % only for corner cases - queue won't grow beyond this. 

%
% Initial conditions: 
%

% 1: initial rate of flow 1
% 2: RTT gradient of flow 1
% ...
% 2*NumFlows: RTT graident of flow numFlowss
% 2*numFlow +1: initial queue size. 
%
initVal = zeros(2*numFlows + 1, 1);

SetInitialRate(1, 5e9);
SetInitialRate(2, 5e9);
SetInitialRate(3, 5e9);

%
% Options.
%
options = ddeset('MaxStep', step_len);

%
% Solve.
%
sol = ddesd(@fluid_timely, @fluid_delays, initVal, [0, sim_length],options);

%
% Extract solution.
%
t = sol.x;
q = sol.y(2*numFlows+1,:);
rates = sol.y(1:2:2*numFlows,:);

%
% Write solution to file.
% 
dlmwrite('timely_fluid.txt',[t',rates'./1e9, q'./8],'\t');

%
% Plot
%
figure
subplot(1,2,1);
plot(t,rates'/1e9);
hold on
axis([0 sim_length 0 10])
xlabel('Time (seconds)')
ylabel('Throughput (Gbps)')

subplot(1,2,2);
plot(t,q./(8e3))
hold on
axis([0 sim_length 0 max(q)/(8e3)])
xlabel('Time (seconds)')
ylabel('Queue (KBytes)')

end

function  SetInitialRate(flownum, rate)
    global initVal;
    initVal(2*(flownum-1)+1, 1) = rate;
end

