function dx = fluid2_test(t,x,lag_matrix)

%{
global target_t;
global last_time;
if t > target_t
    %disp(sprintf('t = %g', target_t))
    target_t = target_t + 0.001
end
last_time = t;
%}

global t0;
global t1;
global Rai;
global C;
global B;
global F;
global g;
global Kmax;
global Kmin;
global timer;
%global timer2;

%ts = d+x(4)/C;

lag = lag_matrix(:,1);

dx = zeros(7,1);

p = h_mark_test(t,lag(7),Kmin,Kmax);
%{
p = 0;
for i=1:2
    if floor(t/t0)==floor((t-t0/2*i)/t0)
        continue
    else
        if h_mark(lag_matrix(7,i),Kmin,Kmax)
           p=1;
        end
    end
end
%}



%p = h_mark(lag(7),Kmin,Kmax);
%p2 = h_mark(lag_matrix(7,2),Kmin,Kmax);
%p3 = h_mark(lag_matrix(7,3),Kmin,Kmax);
%p4 = h_mark(lag_matrix(7,4),Kmin,Kmax);
%p5 = h_mark(lag_matrix(7,5),Kmin,Kmax);

% 7: queue
% 6: alpha2
% 5: rt2
% 4: rc2 
% 3: alpha1
% 2: rt1
% 1: rc1

if p==0
      dx(1)=(x(2)-x(1))*lag(1)/2/B+(x(2)-x(1))/2/timer;
      %if p2==0 && p3==0 && p4==0 && p5==0
      dx(2)=Rai*lag(1)/B+Rai/timer;
      %else
      %    dx(2)=0;
      %end
      dx(3)=-g/t1*x(3);
      dx(4)=(x(5)-x(4))*lag(4)/2/B+(x(5)-x(4))/2/timer;
      %if p2==0 && p3==0 && p4==0 && p5==0
        dx(5)=Rai*lag(4)/B+Rai/timer;
      %else
      %    dx(5)=0;
      %end
      dx(6)=-g/t1*x(6);
else if p==1
            dx(1)=-x(1)*x(3)/2/t0;
            dx(2)=(x(1)-x(2))/t0;
            dx(3)=g/t1*(1-x(3));
            dx(4)=-x(4)*x(6)/2/t0;
            dx(5)=(x(4)-x(5))/t0;
            dx(6)=g/t1*(1-x(6)); 
      else
            a = 1-(1-p)^(t0*lag(1));
            b = p/((1-p)^(-B)-1);
            c = b*((1-p)^(F*B));
            d = p/((1-p)^(-timer*lag(1))-1);
            e = d*((1-p)^(F*timer*lag(1)));
            dx(1) = -x(1)*x(3)/2/t0*a + (x(2)-x(1))*lag(1)/2*b + (x(2)-x(1))*lag(1)/2*d;
            dx(2) = (x(1)-x(2))/t0*a + Rai*lag(1)*c + Rai*lag(4)*e;
            dx(3) = g/t1*((1-(1-p)^(t1*x(1)))-x(3));
            %dx(3) = g/t1*(p-x(3));
            
            a = 1-(1-p)^(t0*lag(4));
            b = p/((1-p)^(-B)-1);
            c = b*((1-p)^(F*B));
            d = p/((1-p)^(-timer*lag(4))-1);
            e = d*((1-p)^(F*timer*lag(4)));
            dx(4) = -x(4)*x(6)/2/t0*a + (x(5)-x(4))*lag(4)/2*b  + (x(5)-x(4))*lag(4)/2*d;
            dx(5) = (x(4)-x(5))/t0*a + Rai*lag(4)*c  + Rai*lag(4)*e;
            dx(6) = g/t1*((1-(1-p)^(t1*x(4)))-x(6));
            %dx(6) = g/t1*(p-x(6));
      end
end

if x(7)>0
    dx(7)=x(4)+x(1)-C;
else
    dx(7)=max(x(4)+x(1)-C,0);
end



