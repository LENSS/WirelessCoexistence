% m=4;
% n=3;
% clear all;
function ret = qos_tuning_unsaturated(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b, phi)

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
% global Wc;
global Epl_b;
global Epl_w;
% global W0;
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
DIFS = 3;

global k;
global scale;
global timeslot;

global Phi;
global Arr_b;
global Arr_w;

global Trw;
global Trb;

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
% W0=aMinBE;%2^aMinBE;
m = 5;%log2(aMaxBE) - log2(aMinBE);
% disp(m);

Trb = ceil((Epl_b)*8*4/30);
Lbt = ceil((Epl_b+phy_b)*8*4/30)+PROP;
Lo = 220;
Wi = Wi_b;
% Wc = Wc_b;

Losw=Los; 
Losb=Lo;    

Epl_b=Epl_b+phy_b;
Epl_w=Epl_w+phy_w;


k = 2;
scale = 3;
timeslot = 1;

Arr_w = arr_w;
Arr_b = arr_b;
Phi = phi;

disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');
disp(['  Phi = ', num2str(Phi) ])
% disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ])
% disp(['  Wc = ', num2str(Wc) ])
disp(['  BMAC packet size = ', num2str(Epl_b) ])
disp(['  BMAC raw packet timeslots = ', num2str(Lbt) ])
disp(['  BMAC arrival rate = ', num2str(Arr_b) ])
% disp(['  aMinBE = ', num2str(2^aMinBE) ])
% disp(['  aMaxBE = ', num2str(2^aMaxBE) ])
disp(['  WiFi packet size = ', num2str(Epl_w) ])
disp(['  WiFi raw packet timeslots = ', num2str(Lsw) ])
disp(['  WiFi arrival rate = ', num2str(Arr_w) ])
disp('======================================================================');
% global z1;

for N_b=BMAC_START:5:BMAC_END
for N_w=WIFI_START:5:WIFI_END
    disp(['  N_w = ', num2str(N_w) ])
    disp(['  N_b = ', num2str(N_b) ])
    disp(['  Total Arr_w = ', num2str(Arr_w*N_w), ' packets/s, or ', num2str(Arr_w*N_w*Epl_w*8/1000000), ' Mbits/s' ]);
    disp(['  Total Arr_b = ', num2str(Arr_b*N_b), ' packets/s, or ', num2str(Arr_b*N_b*Epl_b*8/1000), ' Kbits/s' ]);
% syms Trw;
% syms Lsw;
% syms Lcw;
% syms Trb;
% syms Lbt;
% syms N_w;
% syms N_b;
% 
% Arr_b = 1;
% Arr_w = 10;
% Phi = 2;


for W0_try=[1024 1600 2048]%[64 128 256 512 1024]%5%256%
    for Wc_try=  [30 40 50]%[50 70 90 110 130]%20%[110 130]%
disp('======================================================================');
    
    disp(['  ########Trying W0=', num2str(W0_try)]);
    disp(['  ########Trying Wc=', num2str(Wc_try)]);

y0=[0.4 0.9 0.9 Wc_try W0_try]; 

% Pc = y(1);
% alpha = y(2);
% Pf = y(3);
% Wc = y(4);
% W0 = y(5);


A=[]; %linear inequality constraints (coefficient)
b=[]; %linear inequality constraints (RHS)
Aeq=[]; %linear equality constraints (coefficient)
beq=[]; %linear equality constraints (RHS)

lb = [0,0,0,0,0]; %lower bound 
ub = [1,1,1,8000,4000]; %upper bound 

options = optimset('Algorithm', 'active-set', 'MaxFunEvals', 2000);%, 'FunValCheck', 'on', 'MaxIter', 1000
try
    [y,fval,exitflag] = fmincon(@qos_objfun, y0, A, b, Aeq, beq, lb, ub, @qos_confun, options);
    disp([' exitflag = ', num2str(exitflag)]);
    if exitflag == -2
        disp('  *************No feasible point was found.');
        continue;
    elseif exitflag == 0
        disp('  *************Stopped prematurely.');
        continue;
    end
catch
    disp(['  *************fmincon error!!!']);
    continue;
end


% disp(['y=', num2str(y)]);
% disp(['fval=', num2str(-fval)]);

Pc = y(1);
alpha = y(2);
Pf = y(3);
Wc = y(4);
W0 = y(5);

    disp(['  alpha = ', num2str(alpha) ])
    disp(['  Pc = ', num2str(Pc) ])
    disp(['  Pf = ', num2str(Pf) ])
    disp(['  Wc = ', num2str(Wc) ])
    disp(['  W0 = ', num2str(W0) ])



    disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');

    disp(['  Arr_w = ', num2str(Arr_w*N_w), ' packets/s, or ', num2str(Arr_w*N_w*Epl_w*8/1000000), ' Mbits/s' ]);
    disp(['  Arr_b = ', num2str(Arr_b*N_b), ' packets/s, or ', num2str(Arr_b*N_b*Epl_b*8/1000), ' Kbits/s' ]);


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
    Pb=1-min((Arr_b/100000)/Srv_b,1);
    
    
    disp(['  Pw = ', num2str(Pw) ]);
    disp(['  Pb = ', num2str(Pb) ]);
    
    b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw)*Pc/(1-Pc));
    disp(['  b001 = ', num2str(b001) ]);
%     tau_w = b001/(1-Pc);
    tau_w = Lsw*b001/(1-Pc);

    disp(['  tau_w = ', num2str(tau_w) ]);

    b000 = 1 / (scale*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ scale*Lbt+scale*(Losb-1)+scale/(1-Pb));%
    disp(['  b000 = ', num2str(b000) ]);
%     tau_b = 9*b000;
    tau_b = 3*Lbt*b000;
    disp(['  tau_b = ', num2str(tau_b) ]);

    Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
    Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
    Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
    Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
    Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));

    bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+max(Lcw,scale*Lbt)*Picbw);

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
    
    a = 2*theta_w/delta_w;
    f = @(z) -a.*z.*exp(a.*z);
    Q_w = integral(f,0,200);
    disp(['  Q_w = ', num2str(Q_w) ]);
%     if(Q_w<0) 
%         Q_w=-1;
%     end
    
    delay_w = Q_w/(Arr_w/100000);
    disp(['  delay_w = ', num2str(delay_w*10), ' us'  ]);
    
    
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
    Q_b = integral(f,0,200);
    disp(['  Q_b = ', num2str(Q_b) ]);
%     if(Q_b<0) 
%         Q_b=-1;
%     end
    
    delay_b = Q_b/(Arr_b/100000);
    disp(['  delay_b = ', num2str(delay_b*10), ' us' ]);


    Pi = (1-tau_w)^(N_w)*(1-tau_b)^(N_b);
    Psw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
    Psb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
    Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
    Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
    Pcbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));
    
    
    
    if (Arr_w/100000)/Srv_w > 1
        aggrThr_w(1) = N_w/(mean_w*10^(-5))*Epl_w*8/1000000*(1-Pcbw-Pcw);
    else
            aggrThr_w(1) = Arr_w*N_w*Epl_w*8/1000000;
%         end
    end
    
    if (Arr_b/100000)/Srv_b > 1 %saturated
        aggrThr_b(1) = N_b/(mean_b*10^(-5))*Epl_b*8/1000*(1-Pcbw-Pcb);
    else
        aggrThr_b(1) = Arr_b*N_b*Epl_b*8/1000*(1-Pcbw-Pcb);
    end    
    
   
    disp(['  aggrThr_trps_w = ', num2str(aggrThr_w(1)), ' Mbps']);
    disp(['  aggrThr_trps_b = ', num2str(aggrThr_b(1)), ' Kbps']);
   
    
    Throughput_w = Psw*Trw/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*max(Lcw,scale*Lbt)+Pi);
    Throughput_b = Psb*scale*Trb/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*max(Lcw,scale*Lbt)+Pi);

    Throughput_all = Throughput_w + Throughput_b;

    disp(['  Throughput_all = ', num2str(Throughput_all)]);
    
    
% ===============================old===============================
%     syms tau_w;
%     syms tau_b;
% %     syms tune;
%     
%     tau_w = tau_b/((1-tau_b)/tune+tau_b);
%     
% %     tau_w = 0.017654;
% %     tau_b = 0.0054842;
% %         
% %     r = solve(tau_b/((1-tau_b)/tune+tau_b)-tau_w==0, 'tune');
% %     tune = eval(r);
% %     disp(tune);
% %     
% %     h = tau_b/((1-tau_b)/tune+tau_b)-tau_w;
% %     disp(h);
% %     
% %     return;
% 
% %     Pisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
% %     Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
% %     Picw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
% %     Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
% %     Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));
% %     
% %     bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);
% %     
% %     Throughput_w = bii*(Pisw)*Trw;%Pis/(2-P0);
% %     Throughput_b = 3*Trb*(Pisb)*bii;
%     
%     Pi = (1-tau_w)^(N_w)*(1-tau_b)^(N_b);
%     Psw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
%     Psb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
%     Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
%     Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
%     Pcbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));
%     
%     Ptr_w=1;
%     Psucc_w = Psw*Trw/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*scale*Lbt+Pi);
%     Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
%     
%     
%     Ptr_b=1;
%     Psucc_b = Psb*scale*Trb/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*scale*Lbt+Pi);
%     Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);
% 
%     Throughput_all = Throughput_w + Throughput_b;
% %     disp(Throughput_all);
%     
%     
% %     return;
%     
%     z1 = diff(Throughput_all, tau_b);
%     z2 = simplify(z1);
% %     disp(z2);
%     [num, deno]=numden(z2);
% %     disp(num);
% %     disp(deno);
% % return;
%     zret = solve(num==0, 'tau_b');
% for i=1:length(zret)
%     if isreal(zret(i)) && zret(i)>0 && zret(i)<1
%         disp('  find a tau_b!');
%         tau_b = zret(i);
%         break;
%     end
%     if i==length(zret)
%         error('  no tau_b exist!');
%     end
% end
%     tau_b = eval(tau_b);
%     z2 = eval(z2);
% %     disp(['  z2 = ', num2str(z2) ]);
%     
%     tau_w = eval(tau_w);
% %   tau_w = 0.017654;
% %   tau_b = 0.0054842;
% %   
% %     disp(['  tau_w = ', num2str(tau_w) ]);
% %     disp(['  tau_b = ', num2str(tau_b) ]);
% %     disp(['  N_w = ', num2str(N_w) ])
% %     disp(['  N_b = ', num2str(N_b) ])
% %     disp(['  Pisw = ', num2str(eval(Pisw)) ]);
% %     disp(['  Pisb = ', num2str(eval(Pisb)) ]);
% %     disp(['  Picw = ', num2str(eval(Picw)) ]);
% %     disp(['  Picb = ', num2str(eval(Picb)) ]);
% %     disp(['  Picbw = ', num2str(eval(Picbw)) ]);
% %     disp(['  bii = ', num2str(eval(bii)) ]);
% 
% % return;
% 
% 
%     disp(['  tau_b = ', num2str(tau_b) ]);
%     disp(['  tau_w = ', num2str(tau_w) ]);
% %     Throughput_w = eval(bii)*(eval(Pisw))*Trw;%Pis/(2-P0);
% %     Throughput_b = 3*Trb*(eval(Pisb))*eval(bii);
% %     Psucc_w = Psw*Trw/(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi);
% %     Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
% % 
% %     Psucc_b = Psb*3*Trb/(Psw*Lsw+Psb*3*Lbt+Pcw*Lcw+Pcb*3*Lbt+Pcbw*3*Lbt+Pi);
% %     Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);
% 
%     Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
%     Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
%     Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
%     Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
%     Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));
% 
%     bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+scale*Lbt*Picbw);
% 
%     Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
%     Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
%     Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
%     Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
%     Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));
% 
%     BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+scale*Lbt*Bicb+scale*Lbt*Bisb+scale*Lbt*Bicbw);
% 
%     
%     alpha = 1-BMACii;
%     disp(['  alpha = ', num2str(alpha) ]);
%     Pf = 1-bii;
%     disp(['  Pf = ', num2str(Pf) ]);
%     
% %     Pc = (Lcw*eval(Picw)+3*Lbt*eval(Picbw))*eval(bii);
%     Pc = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b)));
%     disp(['  Pc = ', num2str(Pc) ]);
%     disp(['  Throughput_w = ', num2str(eval(Throughput_w)) ]);
%     disp(['  Throughput_b = ', num2str(eval(Throughput_b)) ]);
%     
% %     Pc = 0.25;
%     
% 
% 
%     timeslot = 1;
% 
%     x = (alpha + (1-alpha)*alpha/k)/1;
% 
%     syms W0;
%     syms Wc;
% 
%     mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
%     mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
%     % mean_b = 3*Lbt - timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(alpha/x - (2*alpha/2*(alpha - 1))/x)*(x - 1) - (timeslot*(3*x - 3)*(Wi - 1)*(x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - (timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(3*x - 3)*(Wc - 1))/2;         
%     % mean_b = 3*Lbt - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - timeslot*(alpha/x - (2*alpha/k*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);
% 
%     Srv_w = 1/(mean_w);
%     Srv_b = 1/(mean_b);
%     %     
%     %     Srv_w = 1/(mean_w)/N_w;
%     %     Srv_b = 1/(mean_b)/N_b;
% 
% %     Pw=1-min((Arr_w/100000)/Srv_w,1);
% %     Pb=1-min((Arr_b/100000)/Srv_b,1);
%     Pw=1-(Arr_w/100000)/Srv_w;
%     Pb=1-(Arr_b/100000)/Srv_b;
% 
% 
%     b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw)*Pc/(1-Pc));
% %     tau_w = Lsw*b001/(1-Pc);
%     % tau_w = b001/(1-Pc);
%     zret1 = solve(Lsw*b001/(1-Pc)==tau_w, 'W0');
%     zret1 = eval(zret1);
% %     disp(['  W0: ', num2str(zret1)]);
%     disp(zret1);
% for i=1:length(zret1)
%     if isreal(zret1(i)) && zret1(i)>0
%         disp('  find a W0!');
%         W0 = zret1(i);
%         break;
%     end
%     if i==length(zret1)
%         error('  no W0 exist!');
%     end
% end
%     disp(['  W0: ', num2str(W0)]);
%     
% 
%     b000 = 1 / (scale*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ scale*Lbt+scale*(Losb-1)+scale/(1-Pb));%
% %     tau_b = Lbt*b000;
%     % tau_b = 9*b000;
%     zret2 = solve(3*Lbt*b000==tau_b, 'Wc');
%     zret2 = eval(zret2);
%     disp(zret2);
% for i=1:length(zret2)
%     if isreal(zret2(i)) && zret2(i)>0
%         disp('  find a Wc!');
%         Wc = zret2(i);
%         break;
%     end
%     if i==length(zret2)
%         error('  no W0 exist!');
%     end
% end
%     disp(['  Wc: ', num2str(Wc)]);
%     %     disp(['  Wc: ', num2str(zret2)]);
% 
%     mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
%     mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
% %     mean_b = 3*Lbt - timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(alpha/x - (2*alpha/2*(alpha - 1))/x)*(x - 1) - (timeslot*(3*x - 3)*(Wi - 1)*(x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - (timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(3*x - 3)*(Wc - 1))/2;         
% %     mean_b = 3*Lbt - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - timeslot*(alpha/x - (2*alpha/k*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);
%     
% %     disp(['  mean_w_all = ', num2str(mean_w*N_w) ]);
% %     disp(['  mean_b_all = ', num2str(mean_b*N_b) ]);
%     disp(['  mean_w = ', num2str(mean_w*10), ' us' ]);
%     disp(['  mean_b = ', num2str(mean_b*10), ' us' ]);
%     disp(['  thr_w = ', num2str(N_w/(mean_w*10^(-5))*Epl_w*8/1000000), ' Mbits/s' ]);
%     disp(['  thr_b = ', num2str(N_b/(mean_b*10^(-5))*Epl_b*8/1000), ' Kbits/s' ]);
%     
%     Srv_w = 1/(mean_w);
%     Srv_b = 1/(mean_b);
% %     
% %     Srv_w = 1/(mean_w)/N_w;
% %     Srv_b = 1/(mean_b)/N_b;
%     
%     disp(['  Srv_w = ', num2str(Srv_w) ]);
%     disp(['  Srv_b = ', num2str(Srv_b) ]);
%     disp(['  Arr_w = ', num2str(Arr_w/100000) ]);
%     disp(['  Arr_b = ', num2str(Arr_b/100000) ]);
% 
%     Lmax = max(Lcw,scale*Lbt);
%     
%     Exp_Xw = - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
%         (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
%     Exp_Fw = Pf/(1-Pf);
%     
%     aw = timeslot+Exp_Fw;
%     bw = Lcw;
%     cw = Exp_Xw;
%     
%    Var_Sw = -(Pc/2 - 1/2)*(Picw*(Lcw + Pf/(Pf - 1))^2 + Picbw*(Lmax + Pf/(Pf - 1))^2 + Pisw*(Lsw + Pf/(Pf - 1))^2 + ...,
%         Picb*(3*Lbt + Pf/(Pf - 1))^2 + Pisb*(3*Lbt + Pf/(Pf - 1))^2)*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
%         (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2*2^m*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))) ...,
%         -(Pc/12 - 1/12)*(timeslot - Pf/(Pf - 1))^2*(1/(Pc - 1) + W0^2/(3*(Pc - 1)) - Pc/(Pc - 1)^2 - (4^(m + 1)*W0^3*(1/(Pc - 1) + ...,
%         (Pc*Pc^m - 1)/(Pc - 1)))/3 + (4*W0^2*(4*4^m*Pc*Pc^m - 1))/(3*(4*Pc - 1))) + ...,
%         aw*cw + (W0*aw^2)/2 + aw^2/4 + cw^2 + (W0^2*aw^2)/4 + ((Pc^2 + Pc)*(aw^2/4 - aw*bw + bw^2))/(Pc - 1)^2 + W0*aw*cw + ...,
%         (Pc*(aw*bw - aw*cw + 2*bw*cw - (W0*aw^2)/2 - aw^2/2 + W0*aw*bw))/(Pc - 1) + ...,
%         ((2*Pc*W0*aw*(aw - 2*bw))/(2*Pc - 1)^2 - (2^(m + 1)*Pc^(m + 1)*W0*aw*(aw - 2*bw)*(m - 2*Pc*m + 1))/(2*Pc - 1)^2)*(Pc - 1) ..., 
%         -(Pc - 1)*((Pc - Pc*Pc^m - Pc*Pc^m*m + Pc^m*Pc^2*m)/(Pc - 1)^2 - Pc/(Pc - 1)^2)*(2^m*W0*aw^2 - 2*2^m*W0*bw*aw) + ...,
%         (Pc - 1)*((W0^2*aw^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - (W0^2*aw^2*(4*4^m*Pc*Pc^m - 1))/(4*Pc - 1) + ...,
%         (W0*aw^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) + (2*W0*aw*cw*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1)) ...,  
%         -(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))*(Pc - 1)*(2^m*W0*aw^2 + 2^m*W0^2*aw^2 - 4^m*W0^2*aw^2 + 2*2^m*W0*aw*cw);
% 
%     disp(['  StdVar_w = ', num2str(sqrt(Var_Sw)) ]);%*N_w^2
%     
%     theta_w = (Arr_w/100000) - Srv_w;
% %     delta_w = (Arr_w/100000)^3*(Arr_w/100000) + (Arr_w/100000)*Srv_w^2*Var_Sw; %Poisson
%     delta_w = 0 + (Arr_w/100000)*Srv_w^2*Var_Sw; %Periodic *N_w^2
%     
%     disp(['  theta_w = ', num2str((theta_w)) ]);
%     disp(['  delta_w = ', num2str(sqrt(delta_w)) ]);
%     
% %     ql_cdf_w(1) = 1- exp(2*theta_w/delta_w*1);
% % %     ql_pdf_w(1) = ql_cdf_w(1);
% % %     disp(['  ql_dist_w = ', num2str((ql_pdf_w(1))) ]);
% %     for queue_length = 2:201
% %         ql_cdf_w(queue_length) = 1- exp(2*theta_w/delta_w*queue_length);
% %         ql_pdf_w(queue_length-1) = ql_cdf_w(queue_length) - ql_cdf_w(queue_length-1);
% % %         disp(['  ql_dist_w = ', num2str((ql_pdf_w(queue_length))) ]);
% %     end
% %     
% % %     plot(ql_pdf_w); hold all;
% %     
% %     mean_ql_w = 0;
% %     for i = 1:200
% %         mean_ql_w = mean_ql_w + ql_pdf_w(i)* i ;
% %     end
% %     
% %     disp(['  mean_ql_w = ', num2str(mean_ql_w) ]);
%     
%     a = 2*theta_w/delta_w;
%     f = @(z) -a.*z.*exp(a.*z);
%     Q_w = integral(f,0,200);
%     disp(['  Q_w = ', num2str(Q_w) ]);
%     if(Q_w<0) 
%         Q_w=-1;
%     end
%     queuelength_w(runNo) = Q_w;
%     
%     delay_w = Q_w/(Arr_w/100000);
%     disp(['  delay_w = ', num2str(delay_w*10), ' us'  ]);
%     totaldelay_w(runNo) = delay_w;
%     
%     
%     Exp_Xb =  - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
% 
%     %     Exp_Xb =  - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - ...,
% %         timeslot*(alpha/x - (2*alpha*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);
%     
%     Exp_Cb = (alpha*timeslot)/x - (2*alpha/k*timeslot*(alpha - 1))/x;
% %     Exp_Cb = (alpha*timeslot)/x - (2*alpha*timeslot*(alpha - 1))/x;
% %     Var_Cb = alpha/x*(timeslot-Exp_Cb)^2+(1-alpha)*alpha/x*(2*timeslot-Exp_Cb)^2;
% 
%     ab = scale*timeslot*(Wi-1)/2;
%     bb = Exp_Cb+scale*timeslot*(Wc-1)/2;
%     cb = Exp_Xb;
% 
%     Var_Sb = (scale*timeslot^2*(Wi^2 - 1))/4 - (x*((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha/k*timeslot*(alpha - 1))/x)^2)/x - ...,
%         (alpha*(alpha - 1)*(2*timeslot - (alpha*timeslot)/x + (2*alpha/k*timeslot*(alpha - 1))/x)^2)/x))/(x - 1) - ...,
%         (scale*timeslot^2*x*(Wc^2 - 1))/(4*(x - 1)) + (ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1);
% 
% 
% 
% 
% %     Var_Sb = - ((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x - (alpha*(alpha - 1)*(2*timeslot - ...,
% %     (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x) - ...,
% %     (3*timeslot^2*(Wc^2 - 1)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/4 - (3*timeslot^2*(Wi^2 - 1)*(x - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/4 + ...,
% %         (ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1);
% 
%     disp(['  StdVar_b = ', num2str(sqrt(Var_Sb)) ]);%*N_b^2
%     
%     
%     theta_b = (Arr_b/100000) - Srv_b;
% %     delta_b = (Arr_b/100000)^3*(Arr_b/100000) + (Arr_b/100000)*Srv_b^2*Var_Sb;%Poisson
%     delta_b = 0 + (Arr_b/100000)*Srv_b^2*Var_Sb;%Periodic *N_b^2
%     
%     disp(['  theta_b = ', num2str((theta_b)) ]);
%     disp(['  delta_b = ', num2str(sqrt(delta_b)) ]);
%     
%     
% %     ql_cdf_b(1) = 1- exp(2*theta_b/delta_b*1);
% % %     ql_pdf_b(1) = ql_cdf_b(1);
% % %     disp(['  ql_dist_w = ', num2str((ql_pdf_b(1))) ]);
% %     for queue_length = 2:201
% %         ql_cdf_b(queue_length) = 1- exp(2*theta_b/delta_b*queue_length);
% %         ql_pdf_b(queue_length-1) = ql_cdf_b(queue_length) - ql_cdf_b(queue_length-1);
% % %         disp(['  ql_dist_b = ', num2str((ql_pdf_b(queue_length))) ]);
% %     end
% %     
% % %     plot(ql_pdf_b); hold all;    
% %     
% %     
% %     
% %     
% %     
% %     mean_ql_b = 0;
% %     for i = 1:200
% %         mean_ql_b = mean_ql_b + ql_pdf_b(i)* i ;
% %     end
% %     
% %     disp(['  mean_ql_b = ', num2str(mean_ql_b) ]);
%     
%     
%     a = 2*theta_b/delta_b;
%     f = @(z) -a.*z.*exp(a.*z);
%     Q_b = integral(f,0,200);
%     disp(['  Q_b = ', num2str(Q_b) ]);
%     if(Q_b<0) 
%         Q_b=-1;
%     end
%     queuelength_b(runNo) = Q_b;
%     
%     delay_b = Q_b/(Arr_b/100000);
%     disp(['  delay_b = ', num2str(delay_b*10), ' us' ]);
%     totaldelay_b(runNo) = delay_b;
% ===============================old===============================


% 
% 
%     syms W0;
%     syms Wc;
%     
%     b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw) +(Lcw+Losw)*Pc/(1-Pc));
% %     b001 = 1 /((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
% %     tau_w = (Lsw)*b001/(1-Pc);
%     
%     zret1 = solve(b001/(1-Pc)==tau_w, 'W0');
% %     zret1 = solve(b001*(1-alpha)/(1-Pc)==tau_w, W0);
%     zret1 = eval(zret1);
%     disp(['  W0: ', num2str(zret1)]);
%     x = alpha + (1-alpha)*alpha;
%     b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*Losb);%
%     
%     zret2 = solve(3*b000==tau_b, 'Wc');
% %     zret1 = solve(b001*(1-alpha)/(1-Pc)==tau_w, W0);
%     zret2 = eval(zret2);
%     disp(['  Wc: ', num2str(zret2)]);
    
    

%     xguess = 0.001;
%     xvect = fsolve('tuning_math', xguess);
%     tau_b = xvect(1);
%     disp(tau_b);
    
    
    
    
% % % %     b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw) +(Lcw+Losw)*Pc/(1-Pc));
% % % %     tau_w = (Lsw)*b001/(1-Pc);
% % % % 
% % % %     x = alpha + (1-alpha)*alpha;
% % % %     b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*Losb);%
% % % %     
% % % %     
% % % % %     Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*(1-P_I^2*(1-(1-phi)^(N_b)));
% % % % %     Pisb = N_b*phi*P_I^2*(1-phi)^(N_b-1)*((1-tau_w)^(N_w-1));
% % % % %     Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*(1-P_I^2*(1-(1-phi)^(N_b)));
% % % % %     Picb = (1 - (1-P_I^2*(1-(1-phi)^(N_b))) - N_b*phi*P_I^2*(1-phi)^(N_b-1))*((1-tau_w)^(N_w-1));
% % % % %     Picbw = (1 - (1-P_I^2*(1-(1-phi)^(N_b))))*(1-((1-tau_w)^(N_w-1)));
% % % % % 
% % % % %     bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);
% % % % 
% % % % %     Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*(1-P_I^2*(1-(1-phi)^(N_b-1)));
% % % % %     Bisb = (N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1)*((1-tau_w)^(N_w));
% % % % %     Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*(1-P_I^2*(1-(1-phi)^(N_b-1)));
% % % % %     Bicb = (1 - (1-P_I^2*(1-(1-phi)^(N_b-1))) - (N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1))*((1-tau_w)^(N_w));
% % % % %     Bicbw = (1 - (1-P_I^2*(1-(1-phi)^(N_b-1))))*(1-((1-tau_w)^(N_w)));
% % % % % 
% % % % %     BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);
% % % %     
% % % %     
% % % %     
% % % %     
% % % %     
% % % %     Ptr_w=1;
% % % %     Psucc_w=bii*(Pisw+Picw+Picbw)*Trw;%Pis/(2-P0);
% % % % %     Psucc_w=bii*(Pisw)*Trw;%Pis/(2-P0);
% % % % 
% % % %     Throughput_w = Ptr_w*Psucc_w;%*Tpl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+ P_idle*Tbo_w);%Ptr_w*Posd_w*Tosd_w);
% % % %     
% % % %     
% % % %     Ptr_b=1;
% % % %     Psucc_b = 3*Trb*(Pisb+Picb+Picbw)*bii;%Pib/(2-P0);
% % % % %     Psucc_b = 3*Trb*(Pisb)*bii;%Pib/(2-P0);
% % % %     Throughput_b = Ptr_b*Psucc_b;%*Tpl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b);
% % % %     
% % % %      aggrThr_w(runNo)=Throughput_w;
% % % %      aggrThr_b(runNo)=Throughput_b;
% % % %     
% % % % 
% % % %     disp(['  aggrThr_w = ', num2str(aggrThr_w(runNo))]);
% % % %     disp(['  aggrThr_b = ', num2str(aggrThr_b(runNo))]);



%     averThr_b(runNo)=;
%     averThr_w(runNo)=;

%     xaxis(runNo)=N;

    runNo = runNo+1;
    end
end
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