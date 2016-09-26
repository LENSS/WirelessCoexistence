% m=4;
% n=3;
% clear all;
function ret = coexist_channel_unsaturated(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b)
% clear all;
% 
% [T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b]=deal(1000000, 20, 20, 40, 40, 16, 1024, 1500, 310, 70, 108);

runNo=1;
% saturated = 0;

% WIFI_START=10;
% WIFI_END=10;
% BMAC_START=15;
% BMAC_END=15;

STEP_SIZE=5;

totalRuns=100;%(WIFI_END-WIFI_START)/STEP_SIZE+1;

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
mactime_w(1:totalRuns)=0;
totaldelay_w(1:totalRuns)=0;
mactime_b(1:totalRuns)=0;
totaldelay_b(1:totalRuns)=0;
P_b(1:totalRuns)=0;
P_w(1:totalRuns)=0;
queuelength_w(1:totalRuns)=0;
queuelength_b(1:totalRuns)=0;

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

global Arr_w;
global Arr_b;

SIFS = 1;
PROP = 1;
ACK = 2;
ACKTO= 2;
DIFS = 3;

global L;
global k;
global scale;

L = 1;%inf;

Epl_b = Packet_b;%20*k_b + 8;
Epl_w = Packet_w;%500*i_w;%1536;%800*i_b-100;%
phy_w=40;%864;%
phy_b=20;
channel_bit_rate = 54;

Trw = ceil((Epl_w)*8/channel_bit_rate/10);
Lsw = ceil((Epl_w+phy_w)*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACK;
Lcw = ceil((Epl_w+phy_w)*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACKTO;


% Ls = 35;
% Lc = 35;
Los = 40;
% W0 = 5;
W0=aMinBE;%2^aMinBE;
m = 5;%log2(aMaxBE) - log2(aMinBE);
% disp(m);

Trb = ceil((Epl_b)*8*4/30);
Lbt = ceil((Epl_b+phy_b)*8*4/30)+PROP;
Lo = 220;
Wi = Wi_b;
Wc = Wc_b;

Losw=Los; 
Losb=Lo;    

Epl_b=Epl_b+phy_b;
Epl_w=Epl_w+phy_w;

Arr_w = arr_w;
Arr_b = arr_b;

% Epl_b=Epl_b+phy_b;
% Epl_w=Epl_w+phy_w;

c_plot = 0; %control plot

%   Wc = 30.6324
%   W0 = 50.7222

% disp(['Lbt: ', num2str(Lbt), ' Lsw: ', num2str(Lsw), ' max: ',num2str(max(Lbt, Lsw))]);
% 
% 
% return;

disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');
disp('***************************unsaturated*********************************');
% disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ]);
disp(['  Wc = ', num2str(Wc) ]);
disp(['  BMAC packet size = ', num2str(Epl_b) ]);
disp(['  BMAC raw packet size = ', num2str(Epl_b) ]);
disp(['  BMAC timeslots per raw packet = ', num2str(Lbt) ]);
disp(['  BMAC arrival rate = ', num2str(Arr_b) ])
disp(['  aMinBE = ', num2str(W0) ]);
disp(['  aMaxBE = ', num2str(aMaxBE) ]);
disp(['  WiFi packet size = ', num2str(Epl_w) ]);
disp(['  WiFi raw packet size = ', num2str(Epl_w) ]);
disp(['  WiFi timeslots per raw packet = ', num2str(Lsw) ]);
disp(['  WiFi arrival rate = ', num2str(Arr_w) ])
disp('======================================================================');

for N_b=BMAC_START:5:BMAC_END
for N_w=WIFI_START:5:WIFI_END
    disp(['  N_w = ', num2str(N_w) ]);
    disp(['  N_b = ', num2str(N_b) ]);
    
    k = 2;
    scale = 3;
    
    
for Arr_b=[1 2 3 4]%0.1:0.1:10%10 %per node arrival rate%5%1%
    
    runNo=1;
    aggrThr_w(1:totalRuns)=0;
    aggrThr_b(1:totalRuns)=0;
    mactime_w(1:totalRuns)=0;
    mactime_b(1:totalRuns)=0;
    totaldelay_w(1:totalRuns)=0;
    totaldelay_b(1:totalRuns)=0;
    P_w(1:totalRuns)=0;
    P_b(1:totalRuns)=0;
    queuelength_w(1:totalRuns)=0;
    queuelength_b(1:totalRuns)=0;
    PCs(1:totalRuns)=0;


for Arr_w=[1 10 20 30]%1:1:100 %per node arrival rate%0.01%30%20%


%     Arr_w = 1/Inter_time_w;%100000000; %number * 10 us
%     Arr_b = 1/Inter_time_b;%10000000;
    
    disp('======================================================================');

    disp(['  Total Arr_w = ', num2str(Arr_w*N_w), ' packets/s, or ', num2str(Arr_w*N_w*Epl_w*8/1000000), ' Mbits/s' ]);
    disp(['  Total Arr_b = ', num2str(Arr_b*N_b), ' packets/s, or ', num2str(Arr_b*N_b*Epl_b*8/1000), ' Kbits/s' ]);

    xguess=[0.1 0.1 0.1]';
    xvect = fsolve('coexist_channel_math_unsaturated', xguess);
    Pc = xvect(1);
    alpha = xvect(2);
    Pf = xvect(3);
    
    disp('The roots from the default "fsolve" are: ')
    disp(['  alpha = ', num2str(alpha) ])
    disp(['  Pc = ', num2str(Pc) ])
    disp(['  Pf = ', num2str(Pf) ])
    PCs(runNo) = Pc;
    
    x = (alpha + (1-alpha)*alpha/k);
%     disp(['  x = ', num2str(x) ])
    timeslot = 1;


    mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
    mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
%     mean_b = 3*Lbt - timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(alpha/x - (2*alpha/2*(alpha - 1))/x)*(x - 1) - (timeslot*(3*x - 3)*(Wi - 1)*(x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - (timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(3*x - 3)*(Wc - 1))/2;         
%     mean_b = 3*Lbt - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - timeslot*(alpha/x - (2*alpha/k*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);
    
%     disp(['  mean_w_all = ', num2str(mean_w*N_w) ]);
%     disp(['  mean_b_all = ', num2str(mean_b*N_b) ]);
    disp(['  mean_w = ', num2str(mean_w*10), ' us' ]);
    disp(['  mean_b = ', num2str(mean_b*10), ' us' ]);
    disp(['  capacity_w = ', num2str(N_w/(mean_w*10^(-5))*Epl_w*8/1000000), ' Mbits/s' ]);
    disp(['  capacity_b = ', num2str(N_b/(mean_b*10^(-5))*Epl_b*8/1000), ' Kbits/s' ]);
    
    Srv_w = 1/(mean_w);
    Srv_b = 1/(mean_b);
%     
%     Srv_w = 1/(mean_w)/N_w;
%     Srv_b = 1/(mean_b)/N_b;
    
    disp(['  Srv_w = ', num2str(Srv_w) ]);
    disp(['  Srv_b = ', num2str(Srv_b) ]);
    disp(['  Arr_w = ', num2str(Arr_w/100000) ]);
    disp(['  Arr_b = ', num2str(Arr_b/100000) ]);
    
    
    Pw=1-min((Arr_w/100000)/Srv_w,1);
    
    if (Arr_w/100000)/Srv_w > 1
        aggrThr_w(runNo) = N_w/(mean_w*10^(-5))*Epl_w*8/1000000;
    else
        aggrThr_w(runNo) = Arr_w*N_w*Epl_w*8/1000000;
    end
    
    
    Pb=1-min((Arr_b/100000)/Srv_b,1);
    
    if (Arr_b/100000)/Srv_b > 1 %saturated
        aggrThr_b(runNo) = N_b/(mean_b*10^(-5))*Epl_b*8/1000;
    else
        aggrThr_b(runNo) = Arr_b*N_b*Epl_b*8/1000;
    end
    
    disp(['  Pw = ', num2str(Pw) ]);
    P_w(runNo) = Pw;
    disp(['  Pb = ', num2str(Pb) ]);
    P_b(runNo) = Pb;
    
    b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw)*Pc/(1-Pc));
    disp(['  b001 = ', num2str(b001) ]);
%     tau_w = b001/(1-Pc);
    tau_w = Lsw*b001/(1-Pc);

    disp(['  tau_w = ', num2str(tau_w) ]);

    b000 = 1 / (scale*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ scale*Lbt+scale*(Losb-1)+scale/(1-Pb));%
    disp(['  b000 = ', num2str(b000) ]);
%     tau_b = 9*b000;
    tau_b = 3*Lbt*b000;
    disp(['  tau_b = ', num2str(tau_b) ]);
    
    
    

    
%=========method old===========
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

    bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+max(Lcw,scale*Lbt)*Picbw);
%     disp(['  bii = ', num2str(bii)]);

    Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
    Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
    Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
    Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
    Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));
    
%     disp(['  Bisw = ', num2str(Bisw) ]);
%     disp(['  Bisb = ', num2str((Bisb)) ]);
%     disp(['  Bicw = ', num2str(Bicw) ]);
%     disp(['  Bicb = ', num2str((Bicb)) ]);
%     disp(['  Bicbw = ', num2str(Bicbw) ]);

    BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+scale*Lbt*Bicb+scale*Lbt*Bisb+max(Lcw,scale*Lbt)*Bicbw);
%     disp(['  BMACii = ', num2str(BMACii)]);
    
    Lmax = max(Lcw,scale*Lbt);
    
    Exp_Xw = - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
        (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
    Exp_Fw = Pf/(1-Pf);
    
    aw = timeslot+Exp_Fw;
    bw = Lcw;
    cw = Exp_Xw;
    
   Var_Sw = -(Pc/2 - 1/2)*(Picw*(Lcw + Pf/(Pf - 1))^2 + Picbw*(Lmax + Pf/(Pf - 1))^2 + Pisw*(Lsw + Pf/(Pf - 1))^2 + ...,
        Picb*(3*Lbt + Pf/(Pf - 1))^2 + Pisb*(3*Lbt + Pf/(Pf - 1))^2)*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
        (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2*2^m*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))) ...,
        -(Pc/12 - 1/12)*(timeslot - Pf/(Pf - 1))^2*(1/(Pc - 1) + W0^2/(3*(Pc - 1)) - Pc/(Pc - 1)^2 - (4^(m + 1)*W0^3*(1/(Pc - 1) + ...,
        (Pc*Pc^m - 1)/(Pc - 1)))/3 + (4*W0^2*(4*4^m*Pc*Pc^m - 1))/(3*(4*Pc - 1))) + ...,
        aw*cw + (W0*aw^2)/2 + aw^2/4 + cw^2 + (W0^2*aw^2)/4 + ((Pc^2 + Pc)*(aw^2/4 - aw*bw + bw^2))/(Pc - 1)^2 + W0*aw*cw + ...,
        (Pc*(aw*bw - aw*cw + 2*bw*cw - (W0*aw^2)/2 - aw^2/2 + W0*aw*bw))/(Pc - 1) + ...,
        ((2*Pc*W0*aw*(aw - 2*bw))/(2*Pc - 1)^2 - (2^(m + 1)*Pc^(m + 1)*W0*aw*(aw - 2*bw)*(m - 2*Pc*m + 1))/(2*Pc - 1)^2)*(Pc - 1) ..., 
        -(Pc - 1)*((Pc - Pc*Pc^m - Pc*Pc^m*m + Pc^m*Pc^2*m)/(Pc - 1)^2 - Pc/(Pc - 1)^2)*(2^m*W0*aw^2 - 2*2^m*W0*bw*aw) + ...,
        (Pc - 1)*((W0^2*aw^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - (W0^2*aw^2*(4*4^m*Pc*Pc^m - 1))/(4*Pc - 1) + ...,
        (W0*aw^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) + (2*W0*aw*cw*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1)) ...,  
        -(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))*(Pc - 1)*(2^m*W0*aw^2 + 2^m*W0^2*aw^2 - 4^m*W0^2*aw^2 + 2*2^m*W0*aw*cw);

    disp(['  StdVar_w = ', num2str(sqrt(Var_Sw)) ]);%*N_w^2
    
    theta_w = (Arr_w/100000) - Srv_w;
%     delta_w = (Arr_w/100000)^3*(Arr_w/100000) + (Arr_w/100000)*Srv_w^2*Var_Sw; %Poisson
    delta_w = 0 + (Arr_w/100000)*Srv_w^2*Var_Sw; %Periodic *N_w^2
    
    disp(['  theta_w = ', num2str((theta_w)) ]);
    disp(['  delta_w = ', num2str(sqrt(delta_w)) ]);
    
%     ql_cdf_w(1) = 1- exp(2*theta_w/delta_w*1);
% %     ql_pdf_w(1) = ql_cdf_w(1);
% %     disp(['  ql_dist_w = ', num2str((ql_pdf_w(1))) ]);
%     for queue_length = 2:201
%         ql_cdf_w(queue_length) = 1- exp(2*theta_w/delta_w*queue_length);
%         ql_pdf_w(queue_length-1) = ql_cdf_w(queue_length) - ql_cdf_w(queue_length-1);
% %         disp(['  ql_dist_w = ', num2str((ql_pdf_w(queue_length))) ]);
%     end
%     
% %     plot(ql_pdf_w); hold all;
%     
%     mean_ql_w = 0;
%     for i = 1:200
%         mean_ql_w = mean_ql_w + ql_pdf_w(i)* i ;
%     end
%     
%     disp(['  mean_ql_w = ', num2str(mean_ql_w) ]);
    


%     1-exp(2*theta_w/delta_w*100000)
    a = 2*theta_w/delta_w;
    f = @(z) -a.*z.*exp(a.*z);
    Q_w = integral(f,0,Inf);
    disp(['  Q_w = ', num2str(Q_w) ]);
%     if(Q_w<0) 
%         Q_w=-1;
%     end
    queuelength_w(runNo) = Q_w;
    
    delay_w = Q_w/(Arr_w/100000);
    disp(['  delay_w = ', num2str(delay_w*10), ' us'  ]);
    totaldelay_w(runNo) = delay_w;
    
    
    Exp_Xb =  - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);

    %     Exp_Xb =  - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - ...,
%         timeslot*(alpha/x - (2*alpha*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);
    
    Exp_Cb = (alpha*timeslot)/x - (2*alpha/k*timeslot*(alpha - 1))/x;
%     Exp_Cb = (alpha*timeslot)/x - (2*alpha*timeslot*(alpha - 1))/x;
%     Var_Cb = alpha/x*(timeslot-Exp_Cb)^2+(1-alpha)*alpha/x*(2*timeslot-Exp_Cb)^2;

    ab = scale*timeslot*(Wi-1)/2;
    bb = Exp_Cb+scale*timeslot*(Wc-1)/2;
    cb = Exp_Xb;

    Var_Sb = (scale*timeslot^2*(Wi^2 - 1))/4 - (x*((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha/k*timeslot*(alpha - 1))/x)^2)/x - ...,
        (alpha*(alpha - 1)*(2*timeslot - (alpha*timeslot)/x + (2*alpha/k*timeslot*(alpha - 1))/x)^2)/x))/(x - 1) - ...,
        (scale*timeslot^2*x*(Wc^2 - 1))/(4*(x - 1)) + (ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1);




%     Var_Sb = - ((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x - (alpha*(alpha - 1)*(2*timeslot - ...,
%     (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x) - ...,
%     (3*timeslot^2*(Wc^2 - 1)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/4 - (3*timeslot^2*(Wi^2 - 1)*(x - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/4 + ...,
%         (ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1);

    disp(['  StdVar_b = ', num2str(sqrt(Var_Sb)) ]);%*N_b^2
    
    
    theta_b = (Arr_b/100000) - Srv_b;
%     delta_b = (Arr_b/100000)^3*(Arr_b/100000) + (Arr_b/100000)*Srv_b^2*Var_Sb;%Poisson
    delta_b = 0 + (Arr_b/100000)*Srv_b^2*Var_Sb;%Periodic *N_b^2
    
    disp(['  theta_b = ', num2str((theta_b)) ]);
    disp(['  delta_b = ', num2str(sqrt(delta_b)) ]);
    
    
%     ql_cdf_b(1) = 1- exp(2*theta_b/delta_b*1);
% %     ql_pdf_b(1) = ql_cdf_b(1);
% %     disp(['  ql_dist_w = ', num2str((ql_pdf_b(1))) ]);
%     for queue_length = 2:201
%         ql_cdf_b(queue_length) = 1- exp(2*theta_b/delta_b*queue_length);
%         ql_pdf_b(queue_length-1) = ql_cdf_b(queue_length) - ql_cdf_b(queue_length-1);
% %         disp(['  ql_dist_b = ', num2str((ql_pdf_b(queue_length))) ]);
%     end
%     
% %     plot(ql_pdf_b); hold all;    
%     
%     
%     
%     
%     
%     mean_ql_b = 0;
%     for i = 1:200
%         mean_ql_b = mean_ql_b + ql_pdf_b(i)* i ;
%     end
%     
%     disp(['  mean_ql_b = ', num2str(mean_ql_b) ]);
    
    
    a = 2*theta_b/delta_b;
    f = @(z) -a.*z.*exp(a.*z);
    Q_b = integral(f,0,Inf);
    disp(['  Q_b = ', num2str(Q_b) ]);
%     if(Q_b<0) 
%         Q_b=-1;
%     end
    queuelength_b(runNo) = Q_b;
    
    delay_b = Q_b/(Arr_b/100000);
    disp(['  delay_b = ', num2str(delay_b*10), ' us' ]);
    totaldelay_b(runNo) = delay_b;
    
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
    
% =================================
    Pi = (1-tau_w)^(N_w)*(1-tau_b)^(N_b);
    Psw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
    Psb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
    Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
    Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
    Pcbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));
    
    
    
    if (Arr_w/100000)/Srv_w > 1
        aggrThr_w(runNo) = N_w/(mean_w*10^(-5))*Epl_w*8/1000000*(1-Pcbw-Pcw);
    else
            aggrThr_w(runNo) = Arr_w*N_w*Epl_w*8/1000000;
%         end
    end
    
    if (Arr_b/100000)/Srv_b > 1 %saturated
        aggrThr_b(runNo) = N_b/(mean_b*10^(-5))*Epl_b*8/1000*(1-Pcbw-Pcb);
    else
        aggrThr_b(runNo) = Arr_b*N_b*Epl_b*8/1000*(1-Pcbw-Pcb);
    end    
    
    
    
    disp(['  Pi = ', num2str((Pi)) ]);
    disp(['  Psw = ', num2str(Psw) ]);
    disp(['  Psb = ', num2str((Psb)) ]);
    disp(['  Pcw = ', num2str(Pcw) ]);
    disp(['  Pcb = ', num2str((Pcb)) ]);
    disp(['  Pcbw = ', num2str(Pcbw) ]);
%     
    
    
    Ptr_w=1;
    Psucc_w=bii*(Pisw+Picw+Picbw)*Lsw;%Pis/(2-P0);
    disp(['  aggrThr_trr_w = ', num2str(Psucc_w)]); %raw data throuput
    Psucc_w=bii*(Pisw+Picw+Picbw)*Trw;%Pis/(2-P0);
    disp(['  aggrThr_trp_w = ', num2str(Psucc_w)]); %payload throuput
%     Psucc_w=bii*(Pisw)*Trw;%Pis/(2-P0);

    Psucc_w = Psw*Trw/(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi);

    Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
    
    
    Ptr_b=1;
    Psucc_b = scale*Lbt*(Bisb+Bicb+Bicbw)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trr_b = ', num2str(Psucc_b)]); %raw data throuput
    Psucc_b = scale*Trb*(Bisb+Bicb+Bicbw)*BMACii;%Pib/(2-P0);
    disp(['  aggrThr_trp_b = ', num2str(Psucc_b)]); %payload throuput
%     Psucc_b = 3*Trb*(Bisb)*BMACii;%Pib/(2-P0);

    Psucc_b = Psb*scale*Trb/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*scale*Lbt+Pi);
%     disp(['  bb = ', num2str(Psb*3*Trb)]); 
%     disp(['  ss = ', num2str(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi)]); 
    
    Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);


%     Throughput_w = Trw*N_w*tau_w*(1-tau_w)^(N_w-1)*(1-tau_b)^N_b;
%     Throughput_b = 3*Trb*N_b*tau_b*(1-tau_b)^(N_b-1)*(1-tau_w)^N_w;
    

%      aggrThr_w(runNo)=Throughput_w;
%      aggrThr_b(runNo)=Throughput_b;
     mactime_w(runNo)= mean_w;
     mactime_b(runNo)= mean_b;
    

    disp(['  aggrThr_trps_w = ', num2str(aggrThr_w(runNo)), ' Mbps']);
    disp(['  aggrThr_trps_b = ', num2str(aggrThr_b(runNo)), ' Kbps']);
    
    
    
end


%     averThr_b(runNo)=;
%     averThr_w(runNo)=;

%     xaxis(runNo)=N;

    runNo = runNo+1;
end

if c_plot
    figure;
    subplot(221);plot(aggrThr_w);xlabel('');ylabel('Througput of WiFi');hold all;
    subplot(222);plot(mactime_w);xlabel('');ylabel('MAC time of WiFi');hold all;
%     subplot(223);plot(totaldelay_w);xlabel('');ylabel('Total delay of WiFi');hold all;
%     subplot(224);plot(P_w);xlabel('');ylabel('Pw of WiFi');hold all;
    subplot(223);plot(queuelength_w);xlabel('');ylabel('Queue length WiFi');hold all;
    subplot(224);plot(PCs);xlabel('');ylabel('Pc');hold all;
    figure;
    subplot(221);plot(aggrThr_b);xlabel('');ylabel('Througput of BMAC');hold all;
    subplot(222);plot(mactime_b);xlabel('');ylabel('MAC time of BMAC');hold all;
%     subplot(224);plot(totaldelay_b);xlabel('');ylabel('Total delay of BMAC');hold all;
    subplot(223);plot(queuelength_b);xlabel('');ylabel('Queue length of BMAC');hold all;
    %subplot(224);plot(P_b);xlabel('');ylabel('Pb of BMAC');hold all;
end
end
end
end




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
