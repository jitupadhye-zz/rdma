function dx = fluid_timely(t,x,lag)
    global callCounter;
    callCounter = callCounter + 1;
    
    dx  = zeros(5, 1);
    
    % 1: rate for flow 1
    % 2: rtt gradiant for flow 1
    % 3: rate for flow 2
    % 4: rtt gradiant for flow 2
    % 5: queue
    
    % lag matrix: 
    % (:,1) is t-t' for flow 1
    % (:,2) is t-t'-t* for flow 1
    % (:,3) is t-t' for flow 2
    % (:,4) is t-t'-t* for flow 2
    
    % update queue delta.
    dx(5) = QueueDelta(x(5), [x(1), x(3)]);
   
    % update rate delta. 
    dx(1) = RateDelta(x(1), lag(5,1), x(2)); 
    dx(3) = RateDelta(x(3), lag(5,3), x(4)); 
    
    % update RTT gradient
    dx(2) = RTTGradientDelta(x(1), x(2), lag(5,1), lag(5,2));
    dx(4) = RTTGradientDelta(x(3), x(4), lag(5,3), lag(5,4));
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
%     if (currentRate < C && deltaRate > 0 && currentRate + deltaRate > C)
%         deltaRate = C - currentRate;
%     end
    
    deltaRate = deltaRate / RTTSampleInterval(currentRate);
    
end

function deltaRTTGradient = RTTGradientDelta(currentRate, currRTTGradient, prevQueue, prevPrevQueue)   
    global alpha;
    global C;
    global minRTT;
    deltaRTTGradient = alpha * (-1 * currRTTGradient + (prevQueue - prevPrevQueue)/(C*minRTT));
    deltaRTTGradient = deltaRTTGradient / RTTSampleInterval(currentRate);
    %deltaRTTGradient = deltaRTTGradient / (Seg/currentRate);
end

function rttSampleInterval = RTTSampleInterval(currentRate)
    global Seg;
    global minRTT;
    %rttSampleInterval = max(Seg/currentRate, minRTT);
    rttSampleInterval = Seg/currentRate;
end





