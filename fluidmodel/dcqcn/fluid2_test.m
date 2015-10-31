function dx = fluid2_test(t,x,lag_matrix)
    global Kmax;
    global Kmin;

    dx = zeros(7,1);
   
    % 1: rc1
    % 2: rt1
    % 3: alpha1
    % 4: rc2
    % 5: rt2
    % 6: alpha2
    % 7: queue
    
    %
    % marking probability
    %
    p = CalculateP(t,lag_matrix(7,1),Kmin,Kmax);
    
    
    %
    % Queue
    %
    dx(7) = QueueDelta(x(7), [x(4), x(1)]);

    %
    % alpha
    %
    dx(3) = AlphaDelta(x(3), lag_matrix(1,1),p);
    dx(6) = AlphaDelta(x(6), lag_matrix(4,1),p);

    %
    % Transmission rate (RC)
    %
    dx(1) = RCDelta(x(1), x(2), x(3), lag_matrix(1, 1), p);
    dx(4) = RCDelta(x(4), x(5), x(6), lag_matrix(4, 1), p);

    %
    % Target rate (RT)
    %
    dx(2) = RTDelta(x(1), x(2), lag_matrix(1, 1), p);
    dx(5) = RTDelta(x(4), x(5), lag_matrix(4, 1), p);
end

% 7: queue
% 6: alpha2
% 5: rt2
% 4: rc2 
% 3: alpha1
% 2: rt1
% 1: rc1

function rcDelta = RCDelta(currRC, currRT, currAlpha, prevRC, p)
    global tau;
    global C;
    %
    % The model cannot correctly handle the case of prevRC = 0. So we stay
    % at 0. 
    % 
    if (currRC == 0 && prevRC == 0)
        rcDelta = 0;
    else 
        [a,b,c,d,e] = IntermediateTerms(p, prevRC);
        rcDelta = -1*currRC*currAlpha*a/(2*tau) + (currRT-currRC)*prevRC*b/2 + (currRT-currRC)*prevRC*d/2;
        %
        % Also, rate cannot exceed C.
        %
        if (prevRC >= C && rcDelta > 0)
            rcDelta = 0;
        end
    end 
end

function rtDelta = RTDelta(currRC, currRT, prevRC, p)
    global tau;
    global Rai;
    % The model cannot correctly handle the case of prevRC = 0. So we stay
    % at 0.
    if (currRC == 0 && prevRC == 0)
        rtDelta = 0;
    else 
        [a,b,c,d,e] = IntermediateTerms(p, prevRC);
        rtDelta = -1*(currRT-currRC)*a/(tau) + Rai*prevRC*c + Rai*prevRC*e;
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

function [a, b, c, d, e] = IntermediateTerms(p, prevRC)
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
                e = d*((1-p)^(F*timer*prevRC));
            end
    end
    if (isnan(a) || isnan(b) || isnan(c) || isnan(d) || isinf(d) || isinf(e))
        fprintf (' *************** ERROR ***********************\n');
        fprintf ('%f ', [a b c d e]);
        fprintf ('\n');
    end
end

