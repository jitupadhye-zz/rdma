function dcqcn2()
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
global numFlows;  % number of flows. 
global initVal; % column of initial values (rate and RTT gradient for each flow, plus initial queue length). 

%
% simulation control.
%
sim_step = 5e-6; % 5 microseconds.
sim_length = sim_step * 20000; 

%
% Fixed parameters. 
%
C = 40 * 1e9;   % 40Gbps. Link speed. 
numFlows = 30;

%
% DCQCN fixed parameters. 
%
tau = 50e-6;  % 50 microseconds. This is the feedback delay. 
tauprime = 55e-6; % 55 microseconds. This is the interval of equation 2. 
F = 5; % Fast recovery steps. 
B = 10 * 8 * 1e6;   %10MB.Byte counter.
Rai = 40 * 1e6; % 40Mbps. Rate increase step.

%
% Tunable parameters.
%
timer = 55e-6; % 55 microseconds. This is rate increase timer.
Kmax = 200 * 8 * 1e3; % 200KB
Kmin = 5 * 8 * 1e3; % 5KB
pmax = 1e-2; % 1 percent. 
g = 1/256; 

%
% Initial conditions: (single column matrix)
%
% 1: rc1
% 2: rt1
% 3: alpha1
% 4: rc2
% 5: rt2
% 6: alpha2
% 7: ...
% 3*numFlows+1: queue
%

initVal = zeros(3*numFlows + 1, 1);

% set traget rate of all flows to C, and alpha to 1
initVal(2:3:3*numFlows,1) = C;
initVal(3:3:3*numFlows,1) = 1;

for i=1:numFlows
    SetInitialRate(i, (1+rand*0.2-0.1)* 1e9);
end

%
% solve 
%
options = ddeset('MaxStep', sim_step);
sol = dde23('fluid2_test', tau, initVal, [0, sim_length], options);

%
% Extract solution.
%
t = sol.x;
q = sol.y(end,:);
rates = sol.y(1:3:end-1,:);

%
% Write solution to file.
% 
dlmwrite('timely_fluid.txt',[t',rates'./1e9, q'./8],'\t');

%
% Plot
%
figure
subplot(3,2,1);
plot(t,rates'/1e9);
hold on
axis([0 sim_length 0 40])
xlabel('Time (seconds)')
ylabel('Throughput (Gbps)')

subplot(3,2,2);
plot(t,q./(8e3))
hold on
axis([0 sim_length 0 max(q)/(8e3)])
xlabel('Time (seconds)')
ylabel('Queue (KBytes)')
end


function  SetInitialRate(flownum, rate)
    global initVal;
    initVal(3*(flownum-1)+1, 1) = rate;
end

