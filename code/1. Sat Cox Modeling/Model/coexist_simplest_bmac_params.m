% m=4;
% n=3;
clear all;

totalRuns = 60;

averThr_w(1:totalRuns)=0;
aggrThr_w(1:totalRuns)=0;
averThr_b(1:totalRuns)=0;
aggrThr_b(1:totalRuns)=0;
alphas(1:totalRuns)=0;
tau_ws(1:totalRuns)=0;
tau_bs(1:totalRuns)=0;

Wis(1:totalRuns)=0;
Wcs(1:totalRuns)=0;
Eplbs(1:totalRuns)=0;
Wis(1:totalRuns)=0;
Wcs(1:totalRuns)=0;
Eplbs(1:totalRuns)=0;
m0s(1:totalRuns)=0;
ms(1:totalRuns)=0;
Eplws(1:totalRuns)=0;

runNo=1;

for k_b=6:6%-1:2  %%%%%Eplb
for j_b=4:4%-1:2  %%%%%Wc
for i_b=4:4%-1:1  %%%%%Wi
for i_w=3:3  %%%%%Eplw
for j_w=4:4  %%%%%m0
for k_w=10:10 %%%%%m+m0
    
    
m0=j_w;
W0=2^m0;
m=k_w-j_w;
n=15;

Wi= 80*i_b;
Wc= 20*j_b;

Losw=60;

Losb=220;

N_b=4;%bmac #
Trials=10;%wifi #
Epl_b = 20*k_b + 8;
Epl_w = 500*i_w+846;%1536;%800*i_b-100;%

Wis(runNo)=Wi;
Wcs(runNo)=Wc;
Eplbs(runNo)=Epl_b;
m0s(runNo)=m0;
ms(runNo)=k_w;
Eplws(runNo)=Epl_w;

disp('***********************************************************************');
disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ])
disp(['  Wc = ', num2str(Wc) ])
disp(['  BMAC packet size = ', num2str(Epl_b) ])
disp(['  aMinBE = ', num2str(m0) ])
disp(['  aMaxBE = ', num2str(k_w) ])
disp(['  WiFi packet size = ', num2str(Epl_w) ])
disp('***********************************************************************');
for N_w=2:Trials

xguess=[0 0 0 Wi Wc Epl_b Epl_w m0 m N_w N_b]';
xvect = fsolve('coexist_math_simplest', xguess);
alpha = xvect(2);
Pc = xvect(1);
Pf = xvect(3);
disp('======================================================================');
disp(['  N_w = ', num2str(N_w) ])
disp(['  N_b = ', num2str(N_b) ])

disp('The roots from the default "fsolve" are: ')
disp(['  alpha = ', num2str(alpha) ])
disp(['  Pc = ', num2str(Pc) ])
disp(['  Pf = ', num2str(Pf) ])

Lbt = ceil(Epl_b*8*4/30);
Lsw = ceil(Epl_w*8/54/10);
Lcw = ceil(Epl_w*8/54/10);
x = alpha + (1-alpha)*alpha;
b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1/ ((W0*(1-((2*Pc)^m+1))/(2*(1-2*Pc)) ) + ( (2^(m+1))*W0*(1-(Pc^(n-m)))*(Pc^(m+1))/(1-Pc) ) + ((1-(Pc^(n+1)))/(2*(1-Pc))) + ((1-(Pc^(n+1)))*(1-Pf)*(Lsw+(Lcw*Pc))/(1-Pc)) );
%tau_w = b001*(1-Pf)*(1-(Pc^(n+1)))/(1-Pc);
tau_w = b001*(1-Pf)/(1-Pc);
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 4*(1-alpha)/(1-x)+3*Lbt+3*Losb);
%b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+3*Lbt+2+3*Losb);
%b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+3*Lbt+2+3*Losb);
tau_b = b000;%/(1-x);

disp(['  tau_w = ', num2str(tau_w) ])
disp(['  tau_b = ', num2str(tau_b) ])

Ptrw = Lsw*(1-Pf)*b001+Lcw*(1-Pf)*Pc*b001/(1-Pc);%(Lsw*b001*(1-Pf) + Lcw*Pc*(1-Pf)*b001/(1-Pc))/(Lsw+Lcw);
Ptrb = 3*Lbt*b000; 
Posw = Losw*b001;
Posb = 3*Losb*b000;

disp(['  Ptrw = ', num2str(Ptrw) ])
disp(['  Ptrb = ', num2str(Ptrb) ])
disp(['  Posw = ', num2str(Posw) ])
disp(['  Posb = ', num2str(Posb) ])

alphas(runNo) = alpha;%(1 - ((1-Ptrw)^N_w)*((1-Ptrb)^(N_b)));
tau_ws(runNo) = tau_w;
tau_bs(runNo) = tau_b;

disp(['  channel busy prob = ', num2str(alphas(runNo)) ])

Ptr_w = 1-((1-tau_w)^N_w);%*(1-tau_b)^N_b;
Psucc_w = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^N_b)/Ptr_w;
Pcolli_w = 1-Psucc_w;
Posd_w = Psucc_w+Pcolli_w;
Pbo_w = (1-Ptrw-Posw)^N_w;

disp(['  Ptr_w = ', num2str(Ptr_w) ])
disp(['  Psucc_w = ', num2str(Psucc_w) ])
disp(['  Pcolli_w = ', num2str(Pcolli_w) ])
disp(['  Posd_w = ', num2str(Posd_w) ])
disp(['  Pbo_w = ', num2str(Pbo_w) ])


Ts_w = Epl_w*8/54;
SIFS = 10;
PROP = 1;
ACK = 40;

DIFS = 60;

Tsucc_w = Ts_w+SIFS+PROP+ACK+DIFS+PROP ;%350;
Tcolli_w = Ts_w+SIFS+PROP+ACK+DIFS+PROP;% 350;
Tbo_w = 10;
Tosd_w = 600;

Tpl_w = Ts_w;

Throughput_w = Ptr_w*Psucc_w*Epl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+Pbo_w*Tbo_w+Ptr_w*Posd_w*Tosd_w);
%Max_tput_w = Ptr_w*Epl_w*8/(Ptr_w*Tsucc_w+Pbo_w*Tbo_w+Ptr_w*Posd_w*Tosd_w);
averThr_w(runNo) = Throughput_w*8/N_w;%/22.6113;%/Max_tput_w;
aggrThr_w(runNo) = Throughput_w*8;%/22.6113;%/Max_tput_w;
% disp(['  Max_throughput_w = ', num2str(Max_tput_w) , ' Mbps'])
disp(['  Throughput_w = ', num2str(Throughput_w) , ' Mbps'])
disp(['  Throughput_w/N = ', num2str(Throughput_w/N_w) , ' Mbps'])
% 
Ptr_b = 1-(1-tau_b)^N_b;%*(1-tau_w)^N_w;
Psucc_b = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^N_w)/Ptr_b;
Pcolli_b = 1-Psucc_b; %???
Posd_b = Psucc_b+Pcolli_b; %???
%Pbo_b = (1-Ptrb-Posb)^N_b;
Pbo_b = (1-Ptrb-Posb)^N_b;

disp(['  Ptr_b = ', num2str(Ptr_b) ])
disp(['  Psucc_b = ', num2str(Psucc_b) ])
disp(['  Pcolli_b = ', num2str(Pcolli_b) ])
disp(['  Posd_b = ', num2str(Posd_b) ])
disp(['  Pbo_b = ', num2str(Pbo_b) ])
% 
Ts_b =  Epl_b*8*4;
Tsucc_b = Ts_b;
Tcolli_b = Ts_b; %???
Tbo_b = 30;
Tosd_b = 220*30;
Tpl_b = Ts_b;

Throughput_b = Ptr_b*Psucc_b*Epl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b+Ptr_b*Posd_b*Tosd_b);
%Max_tput_b = Ptr_b*Epl_b*8000/(Ptr_b*Tsucc_b+Pbo_b*Tbo_b+Ptr_b*Posd_b*Tosd_b);
averThr_b(runNo) = Throughput_b*8000/N_b;%/250;%/Max_tput_b;8*1000
aggrThr_b(runNo) = Throughput_b*8000;%/250;%/Max_tput_b;*8*1000
% disp(['  Max_throughput_b = ', num2str(Max_tput_b) , ' Kbps'])
disp(['  Throughput_b = ', num2str(Throughput_b) , ' Kbps'])
disp(['  Throughput_b/N = ', num2str(Throughput_b/N_b) , ' Kbps'])
disp('======================================================================');
diary on
disp([num2str(Throughput_w), ' ', num2str(Throughput_b)]);
diary off
runNo = runNo+1;
end
end
end
end
end
end
end

plot(aggrThr_w,'DisplayName','aggrThr_w','YDataSource','aggrThr_w');hold all;
plot(averThr_w,'DisplayName','averThr_w','YDataSource','averThr_w');hold all;
plot(aggrThr_b,'DisplayName','aggrThr_b','YDataSource','aggrThr_b');hold all;
plot(averThr_b,'DisplayName','averThr_b','YDataSource','averThr_b');
%plot(alphas,'DisplayName','alpha','YDataSource','alphas');hold all;
%plot(tau_ws,'DisplayName','tau_w','YDataSource','tau_ws');hold all;
%plot(tau_bs,'DisplayName','tau_b','YDataSource','tau_bs');hold all;
hold off;grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,100,0,100]);figure(gcf);
