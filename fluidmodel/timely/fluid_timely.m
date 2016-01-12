function dx = fluid_timely(t,x,lag)
    global numFlows;
    
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
        dx(i) = RateDelta(x(i), lag(2*numFlows+1,i), x(i+1));
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

function deltaRate = RateDelta(currentRate, prevQueue, rttGradient)
    global delta;
    global beta;
    global C;
    global T_high;
    global T_low;
    
    queueLow = C * T_low;
    queueHigh = C * T_high;
    if (prevQueue < queueLow)
       deltaRate = delta;
    else if (prevQueue > queueHigh)
            deltaRate = -1 * beta * (1 - queueHigh/prevQueue) * currentRate;
        else
            if (rttGradient <= 0)
                deltaRate = delta;
            else
                deltaRate = -1 * rttGradient * beta * currentRate;
            end
        end
    end
    
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





