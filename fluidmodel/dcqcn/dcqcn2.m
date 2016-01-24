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
    global packetSize; 
    
    % PI parameters
    global usePI;
    global pold;
    global qold;
    global qref;
    global a;
    global b;
    
    % for debugging purposes.
    global numCalls;
    
    %
    % simulation control.
    %
    sim_step = 500e-6; % 5 microseconds.
    options = ddeset('MaxStep', sim_step, 'RelTol', 1e-2, 'AbsTol', 1e-4);
    sim_length = 500e-3;
    numCalls = 0;
    
    % !!!!!!
    % All quantities (rates, buffers, queues etc.) are specified in units of packet size. The reason
    % is that the probability calculations are per packet.b
    % !!!!
    
    %
    % Fixed parameters.
    %
    packetSize = 4e3;
    C = 40e9/packetSize;   % 40Gbps. Link speed.
    
    % PI parameters
    usePI = 1;
    qold = 0;
    pold = 0;
    qref = 200e3*8/packetSize;
    %    a = 1.822e-3/1.1;
    b = 1.816e-5;
    a = (1+.0001)*b;
    %
    % DCQCN fixed parameters.
    %
    tau = 50e-6;  % 50 microseconds. This is the time 
    %taustar = 4e-6; %1 microsoecond. This is the feedback loop delay.
    tauprime = 55e-6; % 55 microseconds. This is the interval of equation 2.
    F = 5; % Fast recovery steps.
    B = 10e6*8/packetSize;   %10MB.Byte counter.
    Rai = 40e6/packetSize; % 40Mbps. Rate increase step.

    %
    % Tunable parameters.
    %
    timer = 55e-6; % 55 microseconds. This is rate increase timer.
    Kmax = 2000e3*8/ packetSize; % 200KB
    Kmin = 5e3*8/ packetSize; % 5KB
    pmax = 1e-1; % 1 percent.
    g = 1/256;

    for taustar = [85e-6]  % vary the feedback delay
        for numFlows = [64]
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
            sol = dde23(@DCQCNModel, taustar, initVal, [0, sim_length], options);
        
            % Extract solution and write to file.
            t = sol.x;
            q = sol.y(end,:);
            rates = sol.y(1:3:end-1,:);

            [utilization, err] = Utilization(t, rates, q, C);  
            fprintf('utilization: PI=%d flows=%d util=%f err=%d\n', usePI, numFlows, utilization, err);              

            prefix = 'unstable';
            if (usePI == 1)
                prefix = 'unstable.pix';
            end
            fileName = sprintf('%s.%d.%d.dat', prefix, numFlows, taustar*1e6);
            
            % when writing to file, rate is in Gbps, and queue is in KB.
            dlmwrite(fileName,[t',rates'.*packetSize/1e9, q'.*packetSize/8e3], '\t');
            fprintf ('median queue=%f\n', median(q));
            PlotSol(t, q, rates, sim_length, numFlows);
        end
        fclose('all');
    end
end
    
function dx = DCQCNModel(t,x,lag_matrix)
    global Kmax;
    global Kmin;
    global numFlows;
    global numCalls;
    global usePI;
    
    % matrix x:
    % 1: rc1
    % 2: rt1
    % 3: alpha1
    % 4: rc2
    % 5: rt2
    % 6: alpha2
    % ...
    % 3*numFlows+1: queue
    
    % lag matrix: 
    % (:,1) is t-t' for flow 1
    % (:,2) is t-t' for flow 2
    % ....

    dx = zeros(3*numFlows+1,1);
    
    %
    % marking probability
    %
    if (usePI == 1)
        p = CalculatePUsingPI(t, lag_matrix(end,1));
    else
        p = CalculateP(t,lag_matrix(end,1),Kmin,Kmax);
    end
    %
    
    %
    % rates and alpha
    %
    for i = 1:3:3*numFlows
        % The model cannot correctly handle the case of prevRC = 0. So we stay at 0. 
        if (x(i) == 0 && lag_matrix(i, 1) == 0)
            dx(i) = 0;
            dx(i+1) = 0;
        else
            [a,b,c,d,e] = IntermediateTerms(p, lag_matrix(i, 1), x(i), t, i);
            % rc
            dx(i) =  RCDelta(x(i), x(i+1), x(i+2), lag_matrix(i, 1), a, b, d);
            % rt
            dx(i+1) = RTDelta(x(i), x(i+1), lag_matrix(i, 1), a, c, e);
        end
        % alpha
        dx(i+2) = AlphaDelta(x(i+2), lag_matrix(i, 1), p);
    end
    
    %
    % Queue
    %
    rates = x(1:3:3*numFlows);
    dx(end) = QueueDelta(x(end), rates);
    
     numCalls = numCalls +1;
     if (mod(numCalls, 10000) == 0)
         fprintf ('%g %d %f\n', t, numCalls, p);
     end
    
end

function rcDelta = RCDelta(currRC, currRT, currAlpha, prevRC, a, b, d)
    global tau;
    global C;
   
    rcDelta = -1*currRC*currAlpha*a/(2*tau) + (currRT-currRC)*prevRC*b/2 + (currRT-currRC)*prevRC*d/2;
    
    % rate cannot exceed C.
    if (currRC >= C && rcDelta > 0)
        rcDelta = 0;
    end
end

function rtDelta = RTDelta(currRC, currRT, prevRC, a, c, e)
    global tau;
    global Rai;
    global C;
    
    rtDelta = -1*(currRT-currRC)*a/(tau) + Rai*prevRC*c + Rai*prevRC*e;
    
    % rate cannot exceed C.
    if (currRT >= C && rtDelta > 0)
        rtDelta = 0;
    end
    
end

function alphaDelta = AlphaDelta(currentAlpha, prevRC, p)
    global tauprime;
    global g;
    alphaDelta = g * ((1-(1-p)^(tauprime*prevRC))-currentAlpha) / tauprime;
end

function queueDelta = QueueDelta(currentQueue, rates)
    global C;
    queueDelta = sum(rates) - C;
    if (currentQueue <= 0)
        queueDelta = max(0, queueDelta);
    end
    if (currentQueue > 10000*8 && queueDelta > 0)
        queueDelta = 0;
    end
end

function p = CalculateP(t, q, kmin, kmax)
    global pmax;
    if (t >= 0)
        if q <= kmin
            p = 0;
        else if q <= kmax
                p = (q-kmin)/(kmax-kmin)*pmax;
            else % q > lmax
                p = 1;
            end
        end
    else 
        p = 1;
    end 
end

function p = CalculatePUsingPI(t, q)
    global pold;
    global qold;
    global qref;
    global a;
    global b;
    
    %xq = max(q, 0);
    %xqold = max(qold, 0);
    %p = a*(xq - qref) - b*(xqold - qref) + pold;
    %fprintf ('qold=%f pold=%f q=%f p=%f\n', qold, pold, q, p);
    %p = min(max(p, 0), 1);
    %qold  = max(q, 0);
    %pold = p;

    p = a*(q - qref) - b*(qold - qref) + pold;
    pold = p;
    p = min(max(p, 0), 1);
    qold  = q;
end

function [a, b, c, d, e] = IntermediateTerms(p, prevRC, currRC, t, i)
    global tau;
    global B;
    global F;
    global timer;   
    if p == 0
            a = 0;
            b = 1/B;
            c = b;
            if (prevRC == 0)
                d = 0;
            else 
                d = 1/(timer*prevRC);
            end
            e = d;
        else if p == 1
                a = 1;
                b = 0;
                c = 0;
                d = 0;
                e = 0;
            else    
                a = 1-(1-p)^(tau*prevRC);   
                b = p/((1-p)^(-B)-1);
                c = b*((1-p)^(F*B));
                d = p/((1-p)^(-timer*prevRC)-1);
                if (isinf(d))
                    d  = 1/(timer*prevRC);
                    e = d;
                else
                    e = d*((1-p)^(F*timer*prevRC));
                end
            end
    end
    if (isnan(a) || isnan(b) || isnan(c) || isnan(d) || isinf(d) || isinf(e))
        fprintf (' *************** ERROR ***********************\n');
        fprintf ('%i ', [p prevRC currRC t i a b c d e]);
        fprintf ('\n');
        a = 0; 
        b = 1/B;
        if (prevRC == 0)
           d = 0;
        else 
           d = 1/(timer*prevRC);
        end
        e = d;
        
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
    global packetSize;
    figure
    subplot(2,1,1);
    plot(t,rates'*packetSize/1e9);
    hold on
    axis([0 sim_length 0 80/numFlows])
    xlabel('Time (seconds)')
    ylabel('Throughput (Gbps)')
    
    subplot(2,1,2);
    plot(t,q.*packetSize/8e3)
    hold on
    axis([0 sim_length 0 1000])
    xlabel('Time (seconds)')
    ylabel('Queue (KBytes)')
end

% function  SetInitialRate(flownum, rate)
%     global initVal;
%     initVal(3*(flownum-1)+1, 1) = rate;
% end
