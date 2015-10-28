function dx = fluid_timely_old(t,x,lag_matrix)

global Seg;
global delta;
global minRTT;
global beta;
global alpha;
global C;
global T_high;
global T_low;
global callCounter;

callCounter = callCounter + 1;

% t* for the two flows with current rates.
t0(1) = max(Seg/x(1), minRTT);
t0(2) = max(Seg/x(3), minRTT);

lag1 = lag_matrix(:,1); % lag1: t-t' for flow 1
lag2 = lag_matrix(:,2); % lag2: t'-t' - t* for flow 1
lag3 = lag_matrix(:,3); % lag1: t-t' for flow 2
lag4 = lag_matrix(:,4); % lag2: t'-t' - t* for flow 2

dx = zeros(5,1);

%lag matrix:
% 1: rate for flow 1
% 2: rtt gradiant for flow 1
% 3: rate for flow 2
% 4: rtt gradiant for flow 2
% 5: queue

% queue
if x(5)>0
    dx(5)=x(3)+x(1)-C;
else
    dx(5)=max(x(3)+x(1)-C,0);
end

if lag1(5) < C*T_low
    dx(1) = delta/t0(1);
else if lag1(5) > C*T_high
        dx(1) = -beta/t0(1)*(1-C*T_high/lag1(5))*x(1);
    else
        if x(2) <= 0
            dx(1) = delta/t0(1);
        else
            dx(1) = -x(2)*beta/t0(1)*x(1);
        end
    end
end
dx(2) = alpha/t0(1)*(-x(2)+(lag1(5)-lag2(5))/C/minRTT);

if lag3(5) < C*T_low
    dx(3) = delta/t0(2);
else if lag3(5) > C*T_high
        dx(3) = -beta/t0(2)*(1-C*T_high/lag3(5))*x(3);
    else
        if x(4) <= 0
            dx(3) = delta/t0(2);
        else
            dx(3) = -x(4)*beta/t0(2)*x(3);
        end
    end
end
dx(4) = alpha/t0(2)*(-x(4)+(lag3(5)-lag4(5))/C/minRTT);
end