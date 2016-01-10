function sol = dcqcn2()
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
    global taustar;
    global tauprime;
    global numFlows;  % number of flows.
    global initVal; % column of initial values (rate and RTT gradient for each flow, plus initial queue length).
    global numCalls;
    
    %
    % simulation control.
    %
    sim_step = 5e-6; % 5 microseconds.
    options = ddeset('MaxStep', sim_step);
    sim_length = sim_step * 10000;
    numCalls = 0;
    
    %
    % Fixed parameters.
    %
    C = 40 * 1e9;   % 40Gbps. Link speed.

    %
    % DCQCN fixed parameters.
    %
    tau = 50e-6;  % 50 microseconds. This is the time 
    taustar = 2.1e-6; %1 microsoecond. This is the feedback loop delay.
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


    numFlows = 2;
    while (numFlows <= 128)
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

        % set initial and traget rate of all flows to C/N, and alpha to 1
        initVal = zeros(3*numFlows + 1, 1);
        initVal(1:3:3*numFlows,1) = C;
        initVal(2:3:3*numFlows,1) = C;
        initVal(3:3:3*numFlows,1) = 1;

        % solve.
        sol = dde23('fluid2_test', taustar, initVal, [0, sim_length], options);
       
        % Extract solution and write to file.
        t = sol.x;
        q = sol.y(end,:);
        rates = sol.y(1:3:end-1,:);
        [utilization, err] = Utilization(t, rates, q, C);
        
        fprintf('%d %f %d\n', numFlows, utilization, err);       
        fileName =  sprintf('dcqcn.%d.dat', numFlows);
        fileId = fopen (fileName, 'w');
        fprintf(fileId, '## utilization = %f\n', utilization);
        fclose(fileId);
        dlmwrite(fileName,[t',rates'./1e9, q'./8], '-append', 'delimiter','\t');
        
        PlotSol(t, q, rates, sim_length, numFlows);
        numFlows = numFlows * 2;
        break;
    end
end

function [u, err] = Utilization (t, rates, q, C)
    sent = 0;
    tmin = t(1,1);
    tmax = t(1,end);
    max = C * (tmax - tmin);
    err = 0;
    for tindex = 1:(size(t, 2)-1)
        ratesum = 0;
        if (q(tindex) > 1)
            ratesum = C;
        else 
            for flow = 1:size(rates, 1)
                ratesum = ratesum + rates(flow, tindex);
            end
            if (ratesum > C)
                ratesum = C;
                err = err + 1;
            end
        end
        sent = sent + ratesum * (t(1, tindex+1) - t(1, tindex) );
    end
    u = sent/max;
end

function PlotSol(t, q, rates, sim_length, numFlows)
    figure
    subplot(2,1,1);
    plot(t,rates'/1e9);
    hold on
    axis([0 sim_length 0 80/numFlows])
    xlabel('Time (seconds)')
    ylabel('Throughput (Gbps)')
    
    subplot(2,1,2);
    plot(t,q./(8e3))
    hold on
    axis([0 sim_length 0 max(q)/(8e3)])
    xlabel('Time (seconds)')
    ylabel('Queue (KBytes)')
end

% function  SetInitialRate(flownum, rate)
%     global initVal;
%     initVal(3*(flownum-1)+1, 1) = rate;
% end