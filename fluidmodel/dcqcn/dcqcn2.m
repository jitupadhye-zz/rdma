clc;clear all;close all;

sim_length = 0.1;
    
global target_t;
global N;
global t0;
global t1;
global Rai;
global C;
global B;
global F;
global g;
global Kmin;
global Kmax;
global d;
global timer;
global timer2;
global last_time;
global pmax;

t0 = 50/1000000;
t1 = 55/1000000;
Rai = 5000;
C = 5000000;
%B = 150;
B = 300;
F = 5;
g = 1/256;
N = 2;
d = 0.000001;
timer = 1500/1000000;
timer2 = timer;
last_time = 0;

Kmin = 40;
Kmax = 40;
pmax = 1;


kmax_list = [60,180,300,500,1000];
%kmax_list = [40:60:1000];
%kmax_list = [2000];
%pmax_list = [0.1,0.2,0.5,1];

%kmax_list = [60];
pmax_list = [1];


%b_list = [300,10000,100000];
%timer_list = [55,150,1500]./1000000;
b_list = [100000];
timer_list = [55]./1000000;



%timer_list = [55:100:1555]./1000000;

%b_list = 300;
%timer_list = 55./1000000;

% 
% ddeset('AbsTol', 1e-100)
% ddeset('RelTol', 1e-100)

figure

diff=[];

%for i=1:length(kmax_list)
%    for j=1:length(pmax_list)
Kmin = 60;
Kmax = 60;
pmax = 1;
for i=1:length(b_list)
    for j=1:length(timer_list)
       [i,j]
        
        target_t = 0;
        
        B = b_list(i);
        timer = timer_list(j);
%        Kmax = kmax_list(i);
%        pmax = pmax_list(j);
        
        r0 =  C;
        q0 =  0;
        alpha0 = 0;

        options = ddeset('MaxStep', 0.000005);
        %options = [];
        %sol = dde23('fluid2', [d+Kmin/C,2*(d+Kmin/C),3*(d+Kmin/C),4*(d+Kmin/C),5*(d+Kmin/C)], [r0; r0; alpha0; r0/2; r0/2; 0; q0], [0, sim_length],options);
        sol = dde23('fluid2_test', [t0/2], [r0/2; r0; 1; r0/2; r0; 1; q0], [0, sim_length],options);
        %[t0/2]
        %sol = dde23('fluid2', [t0/2,t0*3/2,t0*5/2,t0*7/2,t0*9/2], [r0; r0; alpha0; r0/2; r0/2; 0; q0], [0, sim_length],options);
        %sol = dde23('fluid2', [t0], [r0; r0; 0; r0/2; r0/2; 0; q0], [0, sim_length],options);


        %ddeget(options, 'MaxStep')

        t = sol.x;
        rc = sol.y(1,:);
        rt = sol.y(2,:);
        alpha = sol.y(3,:);
        rc2 = sol.y(4,:);
        rt2 = sol.y(5,:);
        alpha2 = sol.y(6,:);
        q = sol.y(7,:);
        

        dlmwrite(sprintf('fluid%d%d.txt',i,j),[t',rc'./1000000*8,rc2'./1000000*8],'\t');
        
        %subplot(length(pmax_list)*2,length(kmax_list),length(pmax_list)*(i-1)+j);
        subplot(length(b_list),length(timer_list),length(timer_list)*(i-1)+j);
        plot(t,rc./1000000*8,'b');
        hold on
        plot(t,rc2./1000000*8,'r--');
        axis([0 sim_length 0 50])
        set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
        %title(sprintf('Kmax=%d',Kmax));
        title(sprintf('Timer=%.1fus,ByteCounter=%d packets',timer*1000000,B));
        xlabel('Time/Second')
        ylabel('Throughput/Gbps')


%{
        %plot qsize
        subplot(length(pmax_list)*2,length(kmax_list),length(pmax_list)*(i-1)+j+length(pmax_list)*length(kmax_list));
        plot(t,q)
        set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
        title(sprintf('Kmax=%d',Kmax));
        xlabel('Time/Second')
        ylabel('Queue Size/Packets')
        axis([0 sim_length 0 100])
        
        %legend('Flow 1','Flow 2')

        %pause
%}
    end
end

set(gcf,'Units','centimeters','Position',[0 0 35 14])
