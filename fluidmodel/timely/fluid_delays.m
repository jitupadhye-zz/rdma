function delays = fluid_delays(t, x)

global Seg;
global minRTT;
global C;
global prop;
global numFlows;

% x is as follows:
% 1: rate for flow 1
% 2: rtt gradiant for flow 1
% 3: rate for flow 2
% 4: rtt gradiant for flow 2
% ...
% 2*numFlows+1: queue

% delay array is as follows:
% 1: t - t' for flow 1
% 2: t - t' - t* for flow 1
% 3: t - t' for flow 2
% 4: t - t' - t* for flow 2
% ...

delays = zeros(2*numFlows, 1);
tprime = x(end)/C + Seg/C + prop;
for i=1:2:2*numFlows
    tstar = max(Seg/x(i), minRTT);
    delays(i) = t - tprime;
    delays(i+1) = t - tstar - tprime;
end
