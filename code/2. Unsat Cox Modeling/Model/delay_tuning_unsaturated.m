% 
% a = [
% 1;
% 2;
% 4;    
% 3
% ];
% 
% 
% return;

function ret = delay_tuning_unsaturated(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b, phi)
% 
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

% global Arr_w;
% global Arr_b;


N_w = WIFI_START;
N_b = BMAC_START;

% global Pc;
% global Pf;
% global tau_w;

SIFS = 1;
PROP = 1;
ACK = 2;
ACKTO= 2;
DIFS = 3;

global k;
global scale;
global timeslot;


Epl_b = Packet_b;%20*k_b + 8;
Epl_w = Packet_w;%500*i_w;%1536;%800*i_b-100;%
phy_w=40;%846;
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
% aMaxBE = 10;
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
global Arr_b;

Arr_b = arr_b;


% disp(['Lbt: ', num2str(Lbt), ' Lsw: ', num2str(Lsw), ' max: ',num2str(max(Lbt, Lsw))]);
% 
% 
% return;

disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');
% disp(['  run number = ', num2str(runNo) ])
disp(['  Wi = ', num2str(Wi) ]);
% disp(['  Wc = ', num2str(Wc) ]);
disp(['  BMAC packet size = ', num2str(Epl_b) ]);
disp(['  BMAC raw packet timeslots = ', num2str(Lbt) ]);
disp(['  BMAC arrival rate = ', num2str(Arr_b) ])
% disp(['  aMinBE = ', num2str(aMinBE) ]);
% disp(['  aMaxBE = ', num2str(aMaxBE) ]);
disp(['  WiFi packet size = ', num2str(Epl_w) ]);
disp(['  WiFi raw packet timeslots = ', num2str(Lsw) ]);
disp('======================================================================');


    disp(['  N_w = ', num2str(N_w) ])
    disp(['  N_b = ', num2str(N_b) ])



% Arr_b = 1;
for W0_try=[32 64 128 256 512 1024]%5%256%512%
    for Wc_try=[10 30 50 70 90 110 130]%20%[110 130]%
disp('======================================================================');
    
    disp(['  ########Trying W0=', num2str(W0_try)]);
    disp(['  ########Trying Wc=', num2str(Wc_try)]);

y0=[0.4 0.9 0.9 Wc_try W0_try 10]; 


% y0=[0.4 0.9 0.9 10 16 20]; 

% Pc = y(1);
% alpha = y(2);
% Pf = y(3);
% Wc = y(4);
% W0 = y(5);
% Arr_w = y(6);
% Arr_b = y(7);
% 


A=[]; %linear inequality constraints (coefficient)
b=[]; %linear inequality constraints (RHS)
Aeq=[]; %linear equality constraints (coefficient)
beq=[]; %linear equality constraints (RHS)

lb = [0,0,0,0,0,0]; %lower bound 
ub = [1,1,1,2000,1000,50]; %upper bound 

options = optimset('Algorithm', 'active-set');%, 'FunValCheck', 'on', 'MaxFunEvals', 200000, 'MaxIter', 100000
try
    [y,fval,exitflag] = fmincon(@delay_objfun, y0, A, b, Aeq, beq, lb, ub, @delay_confun, options);
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


%  ================================old================================
%   alpha = 0.9;
%   Pc = 0.50174;
%   Pf = 0.96032;
% %   mean_w = 191.2644;
% %   mean_b = 77459.2021;
%   Pw = 0.043678;
%   Pb = 0;
%   tau_w = 0.0058054;
%   tau_b = 3.7722e-05;
%   Wc = 80;
%   W0 = 16;
% Arr_w = 1;%0.005;
% Arr_b = 0.1;%0.000012;
% disp(['  Arr_w = ', num2str(Arr_w) ]);
% disp(['  Arr_b = ', num2str(Arr_b) ]);
%   
% y0=[tau_w,tau_b,alpha,Pb,Wc,Arr_b,Pw,W0,Pf,Pc, Arr_w]; %tau_w tau_b alpha Pb Wc ,Pw,W0,Pf,Pc
% 
% 
% 
% A=[]; %linear inequality constraints (coefficient)
% b=[]; %linear inequality constraints (RHS)
% Aeq=[]; %linear equality constraints (coefficient)
% beq=[]; %linear equality constraints (RHS)
% 
% lb = [0,0,0,0,0,0,0,0,0,0,0]; %lower bound 
% ub = [1,1,1,1,Inf,10,1,Inf,1,1,50]; %upper bound 
% 
% options = optimset('Algorithm', 'active-set', 'MaxFunEvals', 200000, 'MaxIter', 100000);
% [y,fval] = fmincon(@objfun, y0, A, b, Aeq, beq, lb, ub, @confun, options);
% 
% disp(['y=', num2str(y)]);
% disp(['fval=', num2str(-fval)]);
% 
% 
% tau_w = y(1); %tau_w tau_b alpha Pb Wc
% tau_b = y(2);
% alpha = y(3);
% Pb = y(4);
% Wc = y(5);
% Arr_b = y(6);
% Pw = y(6);
% W0 = y(7);
% Pf = y(8);
% Pc = y(9);
% Arr_w = y(10);
%  ================================old================================

Pc = y(1);
alpha = y(2);
Pf = y(3);
Wc = y(4);
W0 = y(5);
Arr_w = y(6);

    disp(['  alpha = ', num2str(alpha) ])
    disp(['  Pc = ', num2str(Pc) ])
    disp(['  Pf = ', num2str(Pf) ])
    disp(['  Wc = ', num2str(Wc) ])
    disp(['  W0 = ', num2str(W0) ])
    disp(['  Arr_w = ', num2str(Arr_w) ])

    disp('======================================================================');

    disp(['  Total Arr_w = ', num2str(Arr_w*N_w), ' packets/s, or ', num2str(Arr_w*N_w*Epl_w*8/1000000), ' Mbits/s' ]);
    disp(['  Total Arr_b = ', num2str(Arr_b*N_b), ' packets/s, or ', num2str(Arr_b*N_b*Epl_b*8/1000), ' Kbits/s' ]);


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
    Q_w = integral(f,0,Inf);
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
    Q_b = integral(f,0,Inf);
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
    disp(['  Throughput_w = ', num2str(Throughput_w)]);

    end
end


return;

% =============================old=============================
% disp(['  tau_w = ', num2str(tau_w) ]);
% disp(['  tau_b = ', num2str(tau_b) ]);
% disp(['  alpha = ', num2str(alpha) ]);
% disp(['  Pb = ', num2str(Pb) ]);
% disp(['  Wc = ', num2str(Wc) ]);
% disp(['  Arr_b = ', num2str(Arr_b) ]);
% 
% 
% Pc = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b)));
% 
% 
%     
% Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
% Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
% Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
% Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
% Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));
% 
% bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+scale*Lbt*Picbw);
% Pf = (1-bii);
%     
% 
% xguess=[0.1 1 1]';
% xvect = fsolve('myfun', xguess);
% Pw = xvect(1);
% W0 = xvect(2);
% Arr_w = xvect(3);
% % 
% % tmp = myfun(xvect);
% % tmp(2)
%     
% disp(['  Arr_w = ', num2str(Arr_w) ]);
% disp(['  Arr_b = ', num2str(Arr_b) ]);
% disp(['  tau_w = ', num2str(tau_w) ])
% disp(['  tau_b = ', num2str(tau_b) ])
% disp(['  alpha = ', num2str(alpha) ])
% disp(['  Pc = ', num2str(Pc) ])
% disp(['  Pf = ', num2str(Pf) ])
% disp(['  Pw = ', num2str(Pw) ]);
% disp(['  Pb = ', num2str(Pb) ]);
% disp(['  Wc = ', num2str(Wc) ]);
% disp(['  W0 = ', num2str(W0) ]);
% 
% 
% 
% 
%     x = (alpha + (1-alpha)*alpha/k);
%     disp(['  x = ', num2str(x) ])
%     timeslot = 1;
% 
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
% 
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
%     
%     Pw=1-min((Arr_w/100000)/Srv_w,1);
%     Pb=1-min((Arr_b/100000)/Srv_b,1);
%     
%     disp(['  Pw = ', num2str(Pw) ]);
%     disp(['  Pb = ', num2str(Pb) ]);
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
%     
%     delay_w = Q_w/(Arr_w/100000);
%     disp(['  delay_w = ', num2str(delay_w*10), ' us'  ]);
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
%     
%     delay_b = Q_b/(Arr_b/100000);
%     disp(['  delay_b = ', num2str(delay_b*10), ' us' ]);
% =============================old=============================


% 
%     disp(['  alpha = ', num2str(alpha) ])
%     disp(['  Pc = ', num2str(Pc) ])
%     disp(['  Pf = ', num2str(Pf) ])
%     
% %     mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
% %     mean_b = 3*Lbt - (timeslot*x*(alpha/x - (2*alpha*(alpha - 1))/x))/(x - 1) + (timeslot*(3*x - 3)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(3*x - 3)*(Wc - 1))/(2*(x - 1)^2);
% %     disp(['  mean_w = ', num2str(mean_w) ]);
% %     disp(['  mean_b = ', num2str(mean_b) ]);
% %     Srv_w = 1/(mean_w);
% %     Srv_b = 1/(mean_b);
% %     Pw=1-Arr_w/Srv_w;
% %     Pb=1-Arr_b/Srv_b;
% %     disp(['  Pw = ', num2str(Pw) ]);
% %     disp(['  Pb = ', num2str(Pb) ]);
%     
%     disp(['  tau_w = ', num2str(tau_w) ]);
% 
%     disp(['  tau_b = ', num2str(tau_b) ]);
% 


% alpha = 0.1;
% Pb=0.2;
% Wc = 100;
% 
% x = sym('alpha + (1-alpha)*alpha');
% tau_b = sym('3 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*(Losb-1)+3/(1-Pb))');%
% %     tau_w = b001/(1-Pc);
% syms tau_w;
% 
% Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
% Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
% Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
% Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
% Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));
% 
% BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);
% % 
% % Pc = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b)));
% % alpha = (1-BMACii);
% % 
%  eval(eval(tau_b));


end