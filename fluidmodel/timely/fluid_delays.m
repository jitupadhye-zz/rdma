function delays = fluid_delays(t, x)

global Seg;
global minRTT;
global C;
global prop;

% 1: rate for flow 1
% 2: rtt gradiant for flow 1
% 3: rate for flow 2
% 4: rtt gradiant for flow 2
% 5: queue

delays = zeros(4,1);

t0 = max(Seg/x(1), minRTT); %t* for flow 1
t1 = max(Seg/x(3), minRTT); %t* for flow 2
t2 = x(5)/C + Seg/C + prop; %t'  for flow 1 and 2

% 1: t - t' for flow 1
% 2: t - t' - t* for flow 1
% 3: t - t' for flow 2
% 4: t - t' - t* for flow 2

delays(1) = t - t2;
delays(2) = t - t0 - t2;
delays(3) = t - t2;
delays(4) = t - t1 - t2;
