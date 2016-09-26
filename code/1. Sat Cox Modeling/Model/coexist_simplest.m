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

WIFI_START=4;
WIFI_END=4;

BMAC_START=1;
BMAC_END=20;

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

Losw=0;%60;

Losb=0;%220;%220;

% N_b=0;%bmac #
% Trials=10;%wifi #
Epl_b = 20*k_b + 8;
Epl_w = 500*i_w;%1536;%800*i_b-100;%
phy_w=846;
phy_b=100;
channel_bit_rate = 54;

Wis(runNo)=Wi;
Wcs(runNo)=Wc;
Eplbs(runNo)=Epl_b;
m0s(runNo)=m0;
ms(runNo)=k_w;
Eplws(runNo)=Epl_w;


disp(' ');
disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');
disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ])
disp(['  Wc = ', num2str(Wc) ])
disp(['  BMAC packet size = ', num2str(Epl_b) ])
disp(['  aMinBE = ', num2str(m0) ])
disp(['  aMaxBE = ', num2str(k_w) ])
disp(['  WiFi packet size = ', num2str(Epl_w) ])
for N_w=WIFI_START:1:WIFI_END
for N_b=BMAC_START:BMAC_END

xguess=[0 0 0 Wi Wc Epl_b+phy_b Epl_w+phy_b m0 m N_w N_b Losw Losb]';
xvect = fsolve('coexist_math_simplest', xguess);
Pc = xvect(1);
alpha = xvect(2);
Pf = xvect(3);
disp('======================================================================');
disp(['  N_w = ', num2str(N_w) ])
disp(['  N_b = ', num2str(N_b) ])

disp('The roots from the default "fsolve" are: ')
disp(['  alpha = ', num2str(alpha) ])
disp(['  Pc = ', num2str(Pc) ])
disp(['  Pf = ', num2str(Pf) ])

Lbt = 1;%ceil((Epl_b+phy_b)*8*4/30);
Lsw = 1;%ceil((Epl_w+phy_w)*8/channel_bit_rate/10);
Lcw = 1;%ceil((Epl_w+phy_w)*8/channel_bit_rate/10);
x = alpha + (1-alpha)*alpha;
b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));

tau_w = b001*(1-Pf)/(1-Pc);
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 4*(1-alpha)/(1-x)+3*Lbt+3*Losb);%

tau_b = b000;%/(1-x);

disp(['  tau_w = ', num2str(tau_w) ])
disp(['  tau_b = ', num2str(tau_b) ])

% Ptrw = Lsw*(1-Pf)*b001+Lcw*(1-Pf)*Pc*b001/(1-Pc);%/(Lsw+Lcw);
% Ptrb = 3*Lbt*b000; 
% Posw = Losw*b001;
% Posb = 3*Losb*b000;
% 
% disp(['  Ptrw = ', num2str(Ptrw) ])
% disp(['  Ptrb = ', num2str(Ptrb) ])
% disp(['  Posw = ', num2str(Posw) ])
% disp(['  Posb = ', num2str(Posb) ])

alphas(runNo) = alpha;
tau_ws(runNo) = tau_w;
tau_bs(runNo) = tau_b;

% disp(['  channel busy prob = ', num2str((1 - ((1-Ptrw)^(N_w))*((1-Ptrb)^(N_b)))) ])

Ptr_w = 1-((1-tau_w)^N_w);%*(1-tau_b)^N_b;
Psucc_w = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^N_b)/Ptr_w;
Pcolli_w = 1-Psucc_w;

%Pbo_w = (1-Ptrw-Posw)^N_w;
%Pos_w = 1-((1-tau_w)^N_w);
P_idle = (1-tau_w)^N_w;% 1-Ptr_w;

disp(['  Ptr_w = ', num2str(Ptr_w) ])
disp(['  Psucc_w = ', num2str(Psucc_w) ])
disp(['  Pcolli_w = ', num2str(Pcolli_w) ])
% disp(['  Posd_w = ', num2str(Pos_w) ])
% disp(['  Pbo_w = ', num2str(Pbo_w) ])
disp(['  P_idle = ', num2str(P_idle) ])

Ts_w = (Epl_w+phy_w)*8/channel_bit_rate;
SIFS = 10;
PROP = 0;
ACK = 20;
%T_ack = (ACK+phy_w)*8/channel_bit_rate;
DIFS = 60;

Tosd_w = 600;
Tsucc_w = Ts_w+SIFS+PROP+ACK+DIFS+PROP+ Tosd_w;%350;
Tcolli_w = Ts_w+SIFS+PROP+ACK+DIFS+PROP+ Tosd_w;% 350;
Tbo_w = 10;

Tpl_w = (Epl_w+phy_w)*8/channel_bit_rate;

Throughput_w = Ptr_w*Psucc_w*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
% Ptr_w = 1-((1-Ptrw)^N_w);%*(1-tau_b)^N_b;
% Psucc_w = N_w*Ptrw*((1-Ptrw)^(N_w-1))*((1-Ptrb)^N_b)/Ptr_w;
% Pcolli_w = 1-Psucc_w;
% Posd_w = Psucc_w+Pcolli_w;
% Pbo_w = (1-Ptrw-Posw)^N_w;
% disp(['  Ptr_w = ', num2str(Ptr_w) ])
% disp(['  Psucc_w = ', num2str(Psucc_w) ])
% disp(['  Pcolli_w = ', num2str(Pcolli_w) ])
% disp(['  Posd_w = ', num2str(Posd_w) ])
% disp(['  Pbo_w = ', num2str(Pbo_w) ])
% Throughput_w = Ptr_w*Psucc_w/(Ptr_w*Psucc_w+Ptr_w*Pcolli_w+Pbo_w+Posw);

Actual_throughput_w = Ptr_w*Psucc_w*(Epl_w+phy_w)*8/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+P_idle*Tbo_w);%Pbo_w*Tbo_w+ Pos_w*Tosd_w);%Ptr_w*Posd_w*Tosd_w);

averThr_w(runNo) = Throughput_w/N_w;
aggrThr_w(runNo) = Throughput_w;
disp(['  Normalized Throughput_w = ', num2str(Throughput_w) , ''])
disp(['  Actual Throughput_w = ', num2str( Actual_throughput_w) , ' Mbps'])
disp(['  Throughput_w/N = ', num2str(Throughput_w/N_w) , ' Mbps'])
% 
Ptr_b = 1-(1-tau_b)^N_b;%*(1-tau_w)^N_w;
Psucc_b = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^N_w)/Ptr_b;
Pcolli_b = 1-Psucc_b; %???
Posd_b = Psucc_b+Pcolli_b; %???
Pbo_b = (1-tau_b)^N_b; %1-Ptr_b;%(1-Ptrb-Posb)^N_b;

disp(['  Ptr_b = ', num2str(Ptr_b) ])
disp(['  Psucc_b = ', num2str(Psucc_b) ])
disp(['  Pcolli_b = ', num2str(Pcolli_b) ])
disp(['  Posd_b = ', num2str(Posd_b) ])
disp(['  Pbo_b = ', num2str(Pbo_b) ])
% 
Tosd_b = 220*30;
Ts_b =  (Epl_b+phy_b)*8*4;
Tsucc_b = Ts_b+Tosd_b;
Tcolli_b = Ts_b+Tosd_b; %???
Tbo_b = 30;

Tpl_b = Ts_b;%(Epl_b)*8*4;

Throughput_b = Ptr_b*Psucc_b*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);

Actual_throughput_b = Ptr_b*Psucc_b*(Epl_b+phy_b)*8000/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);

averThr_b(runNo) = Throughput_b/N_b;%/250;8*1000
aggrThr_b(runNo) = Throughput_b;%/250;*8*1000
disp(['  Normalized Throughput_b = ', num2str(Throughput_b) , ''])
disp(['  Actual Throughput_b = ', num2str(Actual_throughput_b) , ' Kbps'])
disp(['  Throughput_b/N = ', num2str(Throughput_b/N_b) , ' Kbps'])
% diary on
disp([num2str(Throughput_w), ' ', num2str(Throughput_b)]);
% diary off
disp('*******************************end*************************************');
disp(' ');

runNo = runNo+1;
end
end
end
end
end
end
end
end

plot(aggrThr_w,'DisplayName','aggrThr_w','YDataSource','aggrThr_w');hold all;
%plot(averThr_w,'DisplayName','averThr_w','YDataSource','averThr_w');hold all;
plot(aggrThr_b,'DisplayName','aggrThr_b','YDataSource','aggrThr_b');hold all;
%plot(averThr_b,'DisplayName','averThr_b','YDataSource','averThr_b');
%plot(alphas,'DisplayName','alpha','YDataSource','alphas');hold all;
%plot(tau_ws,'DisplayName','tau_w','YDataSource','tau_ws');hold all;
%plot(tau_bs,'DisplayName','tau_b','YDataSource','tau_bs');hold all;
hold off;grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,100,0,1]);figure(gcf);
