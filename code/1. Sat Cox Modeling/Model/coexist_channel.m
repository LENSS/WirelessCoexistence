% m=4;
% n=3;
% clear all;
function ret = coexist_channel(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b)

runNo=1;
% saturated = 0;

% WIFI_START=10;
% WIFI_END=10;
% BMAC_START=15;
% BMAC_END=15;

STEP_SIZE=5;

totalRuns=(WIFI_END-WIFI_START)/STEP_SIZE+1;

xaxis(1:totalRuns)=0;
averThr_w(1:totalRuns)=0;
aggrThr_w(1:totalRuns)=0;
averThr_b(1:totalRuns)=0;
aggrThr_b(1:totalRuns)=0;

alphas(1:totalRuns)=0;
phis(1:totalRuns)=0;
betas(1:totalRuns)=0;
taus(1:totalRuns)=0;
Pfs(1:totalRuns)=0;
PCs(1:totalRuns)=0;

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

% global N;
% global Ls;
% global Lc;
% global Los;
% global W0;
% global m;
% global CWmin;
% global L;
% global tau;

global SIFS;
global PROP;
global ACK;
global ACKTO;
global DIFS;

global Wi;
global Wc;
global Epl_b;
global Epl_w;
global W0;
global m;
global N_w;
global N_b;
global Losw;
global Losb;
global channel_bit_rate;

global Lsw;
global Lcw;
global Lbt;

SIFS = 1;
PROP = 1;
ACK = 2;
ACKTO= 2;
DIFS = 6;

global L;

L = 1;%inf;

Epl_b = Packet_b;%20*k_b + 8;
Epl_w = Packet_w;%500*i_w;%1536;%800*i_b-100;%
phy_w=846;
phy_b=20;
channel_bit_rate = 54;

Trw = ceil((Epl_w)*8/channel_bit_rate/10);
Lsw = ceil((Epl_w+phy_w)*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACK;
Lcw = ceil((Epl_w+phy_w)*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACKTO;


% Ls = 35;
% Lc = 35;
Los = 20;
% W0 = 5;
W0=aMinBE;%2^aMinBE;
m = log2(aMaxBE) - log2(aMinBE);
% disp(m);

Trb = ceil((Epl_b)*8*4/30);
Lbt = ceil((Epl_b+phy_b)*8*4/30)+PROP;
Lo = 220;
Wi = Wi_b;
Wc = Wc_b;

Losw=Los; 
Losb=Lo;    

Epl_b=Epl_b+phy_b;
Epl_w=Epl_w+phy_b;


% disp(['Lbt: ', num2str(Lbt), ' Lsw: ', num2str(Lsw), ' max: ',num2str(max(Lbt, Lsw))]);
% 
% 
% return;

disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');
disp('****************************saturated**********************************');
% disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ])
disp(['  Wc = ', num2str(Wc) ])
disp(['  BMAC packet size = ', num2str(Epl_b) ])
disp(['  BMAC raw packet timeslots = ', num2str(Lbt) ])
disp(['  aMinBE = ', num2str(aMinBE) ])
disp(['  aMaxBE = ', num2str(aMaxBE) ])
disp(['  WiFi packet size = ', num2str(Epl_w) ])
disp(['  WiFi raw packet timeslots = ', num2str(Lsw) ])
disp('======================================================================');

for N_b=BMAC_START:5:BMAC_END
for N_w=WIFI_START:5:WIFI_END
    disp(['  N_w = ', num2str(N_w) ])
    disp(['  N_b = ', num2str(N_b) ])
    
    xguess=[0.1 0.1 0.1]';
    xvect = fsolve('coexist_channel_math', xguess);
    Pc = xvect(1);
    alpha = xvect(2);
    Pf = xvect(3);
    
    disp('The roots from the default "fsolve" are: ')
    disp(['  alpha = ', num2str(alpha) ])
    disp(['  Pc = ', num2str(Pc) ])
    disp(['  Pf = ', num2str(Pf) ])
    

    b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw)+(Lcw+Losw)*Pc/(1-Pc));
    disp(['  b001 = ', num2str(b001) ]);
    tau_w = b001/(1-Pc);
    disp(['  tau_w = ', num2str(tau_w) ]);

    x = alpha + (1-alpha)*alpha;
    b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*Losb);%
    disp(['  b000 = ', num2str(b000) ]);
    
%=========method old===========
    tau_b = 3*b000;
    
    disp(['  tau_b = ', num2str(tau_b) ]);
%     continue;
if N_b == 0
    
    Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1));
    Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)));

    bii = 1/(1+Lsw*Pisw+Lcw*Picw);
    
    Pi = (1-tau_w)^(N_w);
    Psw = N_w*tau_w*((1-tau_w)^(N_w-1));
    Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)));
    
    Ptr_w=1;
    Psucc_w=bii*(Pisw+Picw)*Lsw;%Pis/(2-P0);
    disp(['  aggrThr_trr_w = ', num2str(Psucc_w)]);
    Psucc_w=bii*(Pisw+Picw)*Trw;%Pis/(2-P0);
    disp(['  aggrThr_trp_w = ', num2str(Psucc_w)]);
%     Psucc_w=bii*(Pisw)*Trw;%Pis/(2-P0);

    Psucc_w = Psw*Trw/(Psw*Lsw+Pcw*Lcw+Pi);

    Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);

     aggrThr_w(runNo)=Throughput_w;
%      aggrThr_b(runNo)=Throughput_b;
    

    disp(['  aggrThr_trps_w = ', num2str(aggrThr_w(runNo))]);
%     disp(['  aggrThr_trps_b = ', num2str(aggrThr_b(runNo))]);
    
elseif N_w == 0
    
    Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1));
    Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)));

    BMACii = 1/(1+3*Lbt*Bicb+3*Lbt*Bisb);

    Pi = (1-tau_b)^(N_b);
    Psb = N_b*tau_b*((1-tau_b)^(N_b-1));
    Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)));
    
    Ptr_b=1;
    Psucc_b = 3*Lbt*(Bisb+Bicb)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trr_b = ', num2str(Psucc_b)]);
    Psucc_b = 3*Trb*(Bisb+Bicb)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trp_b = ', num2str(Psucc_b)]);
%     Psucc_b = 3*Trb*(Bisb)*BMACii;%Pib/(2-P0);

    Psucc_b = Psb*3*Trb/(Psb*3*Lbt+Pcb*3*Lbt+Pi);
    Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);


%     Throughput_w = Trw*N_w*tau_w*(1-tau_w)^(N_w-1)*(1-tau_b)^N_b;
%     Throughput_b = 3*Trb*N_b*tau_b*(1-tau_b)^(N_b-1)*(1-tau_w)^N_w;
    

%      aggrThr_w(runNo)=Throughput_w;
     aggrThr_b(runNo)=Throughput_b;
    

%     disp(['  aggrThr_trps_w = ', num2str(aggrThr_w(runNo))]);
    disp(['  aggrThr_trps_b = ', num2str(aggrThr_b(runNo))]);
    
else

    Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
    Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
    Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
    Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
    Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));

    bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);
    disp(['  bii = ', num2str(bii)]);

    Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
    Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
    Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
    Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
    Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));

    BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);
    disp(['  BMACii = ', num2str(BMACii)]);
    
%     disp(['  N_w = ', num2str(N_w) ])
%     disp(['  N_b = ', num2str(N_b) ])
%     disp(['  Pisw = ', num2str(Pisw) ]);
%     disp(['  Pisb = ', num2str(Pisb) ]);
%     disp(['  Picw = ', num2str(Picw) ]);
%     disp(['  Picb = ', num2str(Picb) ]);
%     disp(['  Picbw = ', num2str(Picbw) ]);
%     disp(['  bii = ', num2str(bii) ]);
   
%=========method new===========
%     P_I = 1 - alpha;
%     phi = 3*b000/(1-x);
% 
%     
%     Pisw = (N_w)*tau_w*((1-tau_w)^(N_w-1))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b)));
%     Pisb = 3*Lbt*N_b*phi*P_I^2*(1-phi)^(N_b-1)*((1-tau_w)^(N_w));
%     Picw = (1 - ((1-tau_w)^(N_w)) - (N_w)*tau_w*((1-tau_w)^(N_w-1)))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b)));
%     Picb = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b))) - 3*Lbt*N_b*phi*P_I^2*(1-phi)^(N_b-1))*((1-tau_w)^(N_w));
%     Picbw = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b))))*(1-((1-tau_w)^(N_w)));
% 
%     bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);
    
    
    
%     Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*(1-P_I^2*(1-(1-phi)^(N_b)));
%     Pisb = N_b*phi*P_I^2*(1-phi)^(N_b-1)*((1-tau_w)^(N_w-1));
%     Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*(1-P_I^2*(1-(1-phi)^(N_b)));
%     Picb = (1 - (1-P_I^2*(1-(1-phi)^(N_b))) - N_b*phi*P_I^2*(1-phi)^(N_b-1))*((1-tau_w)^(N_w-1));
%     Picbw = (1 - (1-P_I^2*(1-(1-phi)^(N_b))))*(1-((1-tau_w)^(N_w-1)));
% 
%     bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);

%     Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*(1-P_I^2*(1-(1-phi)^(N_b-1)));
%     Bisb = (N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1)*((1-tau_w)^(N_w));
%     Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*(1-P_I^2*(1-(1-phi)^(N_b-1)));
%     Bicb = (1 - (1-P_I^2*(1-(1-phi)^(N_b-1))) - (N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1))*((1-tau_w)^(N_w));
%     Bicbw = (1 - (1-P_I^2*(1-(1-phi)^(N_b-1))))*(1-((1-tau_w)^(N_w)));
% 
%     BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);
    
%=================================
    Pi = (1-tau_w)^(N_w)*(1-tau_b)^(N_b);
    Psw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
    Psb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
    Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
    Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
    Pcbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));
    
    disp(['  Pi = ', num2str((Pi)) ]);
    disp(['  Psw = ', num2str(Psw) ]);
    disp(['  Psb = ', num2str((Psb)) ]);
    disp(['  Pcw = ', num2str(Pcw) ]);
    disp(['  Pcb = ', num2str((Pcb)) ]);
    disp(['  Pcbw = ', num2str(Pcbw) ]);
    
    Ptr_w=1;
    Psucc_w=bii*(Pisw+Picw+Picbw)*Lsw;%Pis/(2-P0);
    disp(['  aggrThr_trr_w = ', num2str(Psucc_w)]);
    Psucc_w=bii*(Pisw+Picw+Picbw)*Trw;%Pis/(2-P0);
    disp(['  aggrThr_trp_w = ', num2str(Psucc_w)]);
%     Psucc_w=bii*(Pisw)*Trw;%Pis/(2-P0);

    Psucc_w = Psw*Trw/(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi);

    Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
    
    
    Ptr_b=1;
    Psucc_b = 3*Lbt*(Bisb+Bicb+Bicbw)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trr_b = ', num2str(Psucc_b)]);
    Psucc_b = 3*Trb*(Bisb+Bicb+Bicbw)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trp_b = ', num2str(Psucc_b)]);
%     Psucc_b = 3*Trb*(Bisb)*BMACii;%Pib/(2-P0);

    Psucc_b = Psb*3*Trb/(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi);
    Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);


%     Throughput_w = Trw*N_w*tau_w*(1-tau_w)^(N_w-1)*(1-tau_b)^N_b;
%     Throughput_b = 3*Trb*N_b*tau_b*(1-tau_b)^(N_b-1)*(1-tau_w)^N_w;
    

     aggrThr_w(runNo)=Throughput_w;
     aggrThr_b(runNo)=Throughput_b;
    

    disp(['  aggrThr_trps_w = ', num2str(aggrThr_w(runNo))]);
    disp(['  aggrThr_trps_b = ', num2str(aggrThr_b(runNo))]);
end


%     averThr_b(runNo)=;
%     averThr_w(runNo)=;

%     xaxis(runNo)=N;

    runNo = runNo+1;
end
end
% plot(aggrThr_w,'DisplayName','aggrThr_w','YDataSource','aggrThr_w');hold all;
% %plot(averThr_w,'DisplayName','averThr_w','YDataSource','averThr_w');hold all;
% plot(xaxis, aggrThr_b,'DisplayName','aggrThr_b');hold all;
% plot(xaxis, alphas,'DisplayName','alphas');hold all;
% plot(xaxis, betas,'DisplayName','betas');hold all;
% plot(xaxis, phis,'DisplayName','phis');hold all;
% 

% subplot(221);plot(xaxis, taus,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\tau');axis([0,N,0,1]);hold all;%figure(gcf);
% subplot(222);plot(xaxis, PCs,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pc');axis([0,N,0,1]);hold all;%figure(gcf);
% subplot(223);plot(xaxis, Pfs,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pf');axis([0,N,0,1]);%figure(gcf);
% subplot(221);plot(xaxis, aggrtau,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\phi');axis([0,N,0,1]);%figure(gcf);
% subplot(222);plot(xaxis, aggrPc,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('\alpha');axis([0,N,0,1]);%figure(gcf);


% subplot(221);plot(xaxis, phis,'DisplayName','averageCCA','YDataSource','averageCCA');xlabel('Number of nodes');ylabel('\phi');axis([0,N,0,0.07]);%figure(gcf);
% subplot(223);plot(xaxis, alphas,'DisplayName','averageBusy1','YDataSource','averageBusy1');xlabel('Number of nodes');ylabel('\alpha');axis([0,N,0,1.3]);%figure(gcf);
% subplot(222);plot(xaxis, betas,'DisplayName','averageBusy2','YDataSource','averageBusy2');xlabel('Number of nodes');ylabel('\beta');axis([0,N,0,1]);%figure(gcf);
% subplot(224);plot(xaxis, taus,'DisplayName','aggregateThr','YDataSource','aggregateThr');xlabel('Number of nodes');ylabel('\tau');axis([0,50,0,0.02]);%figure(gcf);
% subplot(224);plot(xaxis, aggrThr_b,'DisplayName','aggregateThr','YDataSource','aggregateThr');xlabel('Number of nodes');ylabel('Thr');%axis([0,50,0,1]);%figure(gcf);

% %plot(averThr_b,'DisplayName','averThr_b','YDataSource','averThr_b');
% %plot(alphas,'DisplayName','alpha','YDataSource','alphas');hold all;
% %plot(tau_ws,'DisplayName','tau_w','YDataSource','tau_ws');hold all;
% %plot(tau_bs,'DisplayName','tau_b','YDataSource','tau_bs');hold all;
% hold off;grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,100,0,1]);figure(gcf);
end