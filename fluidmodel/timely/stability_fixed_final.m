clear all; close all; clc

N = 2;
C = 10 * 1e9; % line rate.
Seg = 16 * 8 * 1e3; % burstsize.
MTU = 1 * 8 * 1e3;
prop = 100 /1e6; % propagation delay

%
% Parameters we can play with.
%
delta = 10 * 1e6; % 10Mbps
T_high = 500/1e6; % 500 microseconds (see section 4.4)
T_low = 50/1e6; % 50 microseconds (see section 4.4). 
minRTT = 20/1e6; % 20 microseconds 
beta = 8/1000;
alpha = 0.875; % unsure
ws = 0.5;
wk = 0.5;

%Rai = 20 * 1e6;
%N = 24;
%T = timer * 2;
%B = B * 10;
%sweep = 1e6 * [30:10:100];
sweep = 2:2:50;
%sweep = [2:2:20]./1000;  %beta
%sweep = [2:2:30] .* 1e6; %delta
%sweep = [0.1:0.1:1]; %ws or alpha
%sweep = [0.1:0.1:1]; %scale
%sweep = [8:8:128] * 8 * 1e3; % Seg
pdiff = zeros(1,length(sweep));

for i=1:length(sweep)
    
    N = sweep(i)
    
    qref = C*T_low;
    qs = N*delta*qref/beta/C*(1-ws)/ws + qref;
    Rs = C / N;
    taus = Seg/Rs;
    taup = qs / C + MTU/C + prop;
    
    s = tf('s');
    
    numerator1 = 2*Rs*delta*alpha/s/minRTT/(Seg*s+alpha*Rs)*wk+Rs*Rs*ws*beta/qref*N/s;
    numerator2 = -2*Rs*delta*alpha/s/minRTT/(Seg*s+alpha*Rs)*wk;    
    denominator = s*Seg+(3*ws-1)*delta;
    hd = numerator1 / denominator * exp(-s*taup) + numerator2 / denominator * exp(-s*(taus+taup));
    [Gm,Pm,Wgm,Wpm] = margin(hd);
    
    pdiff(i) = Pm;

    
end

figure
plot(sweep,pdiff)
hold on
%plot(sweep,pdiff2, 'r')
%plot(sweep,pdiff3, 'g')
ylabel('Phase margin')
xlabel('N')
%legend('2','8','16')

dlmwrite('pm_timely_delay100_fixed.txt',[sweep',pdiff'], 'delimiter','\t');


