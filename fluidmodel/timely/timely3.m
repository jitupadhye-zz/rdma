function sol = timely3()
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
    global sim_length;
    global a; % PI parameters
    global b; % PI parameters
    global prevprevQueue %hack
    global pold;
    global p;
    global qold;
    pold=0;
    p=0;
    prevprevQueue = 0;
    qold = 0;

    %
    % Simulation control
    % 
    step_len = 5e-6 ; % 5 microseconds
    sim_length = 5e-2; % 50 milliseconds 

    % 
    % Fixed Parameters
    %
    C = 10 * 1e9; % line rate.
    Seg = 64 * 8 * 1e3; % burstsize.
    prop = 4e-6; % propagation delay
    
    % Setting PI parameters
    a = (1.822e-9);
    b = (1.82199998e-9);

    %
    % Parameters we can play with.
    %
    delta = 10e8; % 10Mbps
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

    for numFlows = 3
        initVal = zeros(2*numFlows + 1, 1);
        for i=1:numFlows
            %SetInitialRate(i, (1+rand*0.2-0.1)* 1e9);
            %SetInitialRate(i, C/numFlows);
        end

        SetInitialRate(1, 1e9);
        SetInitialRate(2, 1e9);
        SetInitialRate(3, 1e9);
        %
        % Options.
        %
        options = ddeset('MaxStep', step_len, 'RelTol', 1e-2, 'AbsTol', 1e-5);

        %
        % Solve.
        %
        sol = ddesd(@TimelyModel, @DelayModel, initVal, [0, sim_length],options);

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
        dlmwrite(fileName,[t',rates'./1e9, q'./8e3], '-append', 'delimiter','\t');
  
       
        PlotSol(t, q, rates);
        %break;
    end
end

function dx = TimelyModel(t,x,lag)
    global numFlows;
    global prevprevQueue;
    dx  = zeros(2*numFlows+1, 1);
    
    % 1: rate for flow 1
    % 2: rtt gradiant for flow 1
    % 3: rate for flow 2
    % 4: rtt gradiant for flow 2
    % ...
    % 2*numFlows+1: queue
    
    % lag matrix: 
    % (:,1) is t-t' for flow 1
    % (:,2) is t-t'-t* for flow 1
    % (:,3) is t-t' for flow 2
    % (:,4) is t-t'-t* for flow 2
    % ....
     
    rates = x(1:2:2*numFlows);
    dx(end) = QueueDelta(x(end), rates);
    % update rate delta. 
    for i = 1:2:2*numFlows
        dx(i) = RateDelta(x(i), lag(2*numFlows+1,i), x(i+1), ...
                          prevprevQueue);
        prevprevQueue = lag(2*numFlows+1,i);
    end  
    
    % update RTT gradient
    for i = 2:2:2*numFlows
        dx(i) = RTTGradientDelta(x(i-1), x(i), lag(2*numFlows+1,i-1), lag(2*numFlows+1,i));
    end  
end

function deltaQueue = QueueDelta(currentQueue, flowRates)
    global C;
    global maxQueue;
    if (currentQueue > 0)
        if (currentQueue < maxQueue)
            deltaQueue = sum(flowRates)-C;
        else 
            deltaQueue = min(sum(flowRates)-C, 0);
        end
    else
        deltaQueue = max(sum(flowRates)-C, 0);
    end
end

function deltaRate = RateDelta(currentRate, prevQueue, rttGradient, ...
                               prevprevQueue)
    global delta;
    global beta;
    global C;
    global T_high;
    global T_low;
    global a;
    global b;
    global pold;
    global p;
    global qold;
    queueLow = C * T_low;
    %queueLow = 0;
    queueHigh = C * T_high;
    qref = (queueHigh+queueLow)/2;
    error = prevQueue-prevprevQueue;

    %p = a*(prevQueue - qref) - b*(prevprevQueue - qref) + pold;
    %p = min(max(p, 0), 1);
    %qold  = prevQueue;
    %pold = p;
    %deltaRate = delta - p*currentRate;

    p = a*(prevQueue - qref) - b*(qold - qref) + pold
    p = min(max(p, 0), 1);
    qold  = prevQueue;
    pold = p;
    deltaRate = delta - p*beta*currentRate
    
    %   if (prevQueue < queueLow)
        %deltaRate = delta;
        %    deltaRate = -1 * gradient * b* currentRate ...
        %   +1*error*currentRate*a;
        % deltaRate
        %else if (prevQueue > queueHigh)
            % deltaRate = -1 * beta * (1 - queueHigh/prevQueue) *
            % currentRate;
            %        deltaRate = -1 * gradient * b* currentRate ...
            %       +1*error*currentRate*a;
            %else
            %if (rttGradient < 0)
                %deltaRate = delta;
                %   deltaRate = -1 * gradient * b * currentRate ...
                %   +1*error*currentRate*a;
                %deltaRate
                %else
                %deltaRate = -1 * gradient * b* currentRate ...
                %   +1*error*currentRate*a;
                %deltaRate
                %deltaRate = -1 * rttGradient * beta * currentRate;
                %end
                %end
                %end
    
    % do not exceed line rate.
    if (currentRate >= C && deltaRate > 0)
        deltaRate = 0;
    end
    deltaRate = deltaRate / RTTSampleInterval(currentRate);
end

function deltaRTTGradient = RTTGradientDelta(currentRate, currRTTGradient, prevQueue, prevPrevQueue)   
    global alpha;
    global C;
    global minRTT;
    deltaRTTGradient = alpha * (-1 * currRTTGradient + (prevQueue - prevPrevQueue)/(C*minRTT));
    deltaRTTGradient = deltaRTTGradient / RTTSampleInterval(currentRate);
end

function rttSampleInterval = RTTSampleInterval(currentRate)
    global Seg;
    rttSampleInterval = Seg/currentRate;
end

function delays = DelayModel(t, x)
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

function PlotSol(t, q, rates)
    global C;
    global sim_length;
    
    figure
    subplot(2,1,1);
    plot(t,rates'/1e9);
    hold on
    axis([0 sim_length 0 C/1e9])
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

