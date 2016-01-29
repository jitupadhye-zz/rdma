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
global Rcs;
    
N = 2;
C = 100 * 1e9 / 8e3;   % 40Gbps. Link speed. 
tau = 50 * 1e-6;  % 50 microseconds. This is the feedback delay. 
tauprime = 55 * 1e-6; % 55 microseconds. This is the interval of equation 2. 
F = 5; % Fast recovery steps. 
B = 10 * 8 * 1e6  / 8e3;   %10MB.Byte counter.
Rai = 10 * 1e6  / 8e3; % 40Mbps. Rate increase step.
timer = 55 * 1e-6; % 55 microseconds. This is rate increase timer.
T = timer;
Kmax = 1000 * 8 * 1e3  / 8e3; % 200KB
Kmin = 5 * 8 * 1e3  / 8e3; % 5KB
pmax = 1e-2; % 1 percent. 
g = 1/256; 
%taus = tau / 2;
taus = 100 * 1e-6;

%Rai = 20 * 1e6;
%N = 24;
%T = timer * 2;
%B = B * 10;
%sweep = 1e6 * [30:10:100];
sweep = 2:2:100;

pdiff = zeros(1,length(sweep));
global Rcs;

for i=1:length(sweep)
    
    %Rai = sweep(i);
    N = sweep(i)
    Rcs = C / N;
    
    %ps = fzero('solvep', 0.0001);
    
    ps = (Rai*N*N/tauprime/C/C*((1/B+N/C/T)^2))^(1/3)
    
    %x = 1-(1-ps)^(Rcs*tauprime)
        
    % ps
    as = 1 - (1-ps)^(tau*Rcs);
    cs = (1-ps)^(F*B)*ps / ((1-ps)^(-B)-1);
    es = (1-ps)^(F*T*Rcs)*ps / ((1-ps)^(-T*Rcs)-1);
    Rts = C / N * (1 + (cs+es)*tau*Rai/as);
    qs = ps/pmax*(Kmax-Kmin) + Kmin;
    alphas = 1 - (1-ps)^(tauprime*C/N);
    A = (1/B + 1/T/Rcs);

    % e^(-s*taus)/s --> e^(-s*taus)*(Ns+gC)*(s+ps*Rcs)
    AA = (-1/2*Rcs*Rcs*alphas-(1/2+A/4)*Rcs*Rts+(1/2+A/4)*Rcs*Rcs)*N*pmax/(Kmax-Kmin);
    % e^(-s*taus) --> e^(-s*taus)*(Ns+gC)*(s+ps*Rcs)*s
    BB = -1/2*ps*Rcs*alphas-A/2*Rcs+(1/2+A/4)*ps*Rcs;
    % 1 --> (Ns+gC)*(s+ps*Rcs)*s
    CC = -1/2*ps*Rcs*alphas-A/2*Rcs+A/2*Rts+(1/2+A/4)*ps*Rcs-(1/2+A/4)*ps*Rts;
    % e^(-s*taus)/s/(Ns+gC) --> e^(-s*taus)*(s+ps*Rcs)
    DD = -1/2*ps*Rcs*Rcs*g*C*N*pmax/(Kmax-Kmin);
    % 1/(s+ps*Rcs) --> (Ns+gC)*s
    EE = (A/2*Rcs-(1/2+A/4)*ps*Rcs)*(ps*Rcs+A*Rai);
    % e^(-s*taus)/(s+ps*Rcs) --> e^(-s*taus)*(Ns+gC)*s
    FF = (A/2*Rcs-(1/2+A/4)*ps*Rcs)*(-ps*Rts+ps*Rcs+A*Rai-(1+2*F+A/2)*Rai*ps);
    % e^(-s*taus)/(s+ps*Rcs)/s --> e^(-s*taus)*(Ns+gC)
    GG = (A/2*Rcs-(1/2+A/4)*ps*Rcs)*(-Rcs*Rts+Rcs*Rcs-(1+2*F+A/2)*Rai*Rcs)*N*pmax/(Kmax-Kmin);

    %syms w;
    %s = 1i*w;
    %numerator = expand((N*s+g*C)*(s+ps*Rcs)*AA + (N*s+g*C)*(s+ps*Rcs)*s*BB + (s+ps*Rcs)*DD + (N*s+g*C)*s*FF + (N*s+g*C)*GG);
    %denominator = expand(-s*(N*s+g*C)*(s+ps*Rcs)*s+(N*s+g*C)*(s+ps*Rcs)*s*CC+s*(N*s+g*C)*EE);
    %wgc = vpa(solve(abs(numerator/denominator)-1, w));
    %wpc = vpa(solve(angle(numerator/denominator*exp(-1i*w*taus))-pi, w));
    %pdiff(i) = wgc(1) - wpc(1);
    
    s = tf('s');
    numerator = ((N*s+g*C)*(s+ps*Rcs)*AA + (N*s+g*C)*(s+ps*Rcs)*s*BB + (s+ps*Rcs)*DD + (N*s+g*C)*s*FF + (N*s+g*C)*GG);
    denominator = (-s*(N*s+g*C)*(s+ps*Rcs)*s+(N*s+g*C)*(s+ps*Rcs)*s*CC+s*(N*s+g*C)*EE);
    hd = numerator / denominator * exp(-s*taus);
    [Gm,Pm,Wgm,Wpm] = margin(hd);
    pdiff(i) = Pm;
end

dlmwrite('pm_100gbps_delay100_fixed.txt',[sweep',pdiff'], 'delimiter','\t');

plot(sweep,pdiff)
ylabel('Phase margin')
xlabel('N')

