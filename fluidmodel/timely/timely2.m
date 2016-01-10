function sol = timely2()
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
    global totrate;
    
    %
    % Simulation control
    % 
    step_len = 5e-6 ; % 5 microseconds
    sim_length = step_len * 10000 ;

    % 
    % Fixed Parameters
    %
    C = 10 * 1e9; % line rate.
    Seg = 64 * 8 * 1e3; % burstsize.
    prop = 4e-6; % propagation delay
    

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

    numFlows = 2;
    while (numFlows < 128)
   
        initVal = zeros(2*numFlows + 1, 1);
        for i=1:numFlows
            %SetInitialRate(i, (1+rand*0.2-0.1)* 1e9);
            SetInitialRate(i, C/numFlows);
        end

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
        [utilization, err] = Utilization(t, rates, q, C);
        
        %
        % Write solution to file.
        % 
        
        fprintf('%d %f %d\n', numFlows, utilization, err);       
        fileName =  sprintf('timely.%d.dat', numFlows);
        fileId = fopen (fileName, 'w');
        fprintf(fileId, '## utilization = %f\n', utilization);
        fclose(fileId);
        dlmwrite(fileName,[t',rates'./1e9, q'./8], '-append', 'delimiter','\t');
  
        numFlows = numFlows * 2; 
        
        %PlotSol(t, q, rates);
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

function PlotSol(t, q, rates, sim_length)
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

function  SetInitialRate(flownum, rate)
    global initVal;
    initVal(2*(flownum-1)+1, 1) = rate;
end

