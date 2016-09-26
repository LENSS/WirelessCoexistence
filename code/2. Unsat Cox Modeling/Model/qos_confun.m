function [c,ceq] = qos_confun(y) %tau_w and Wc

global Wi;
% global Wc;
% global Epl_b;
% global Epl_w;
% global W0;
global m;
global N_w;
global N_b;
global Losw;
global Losb;
% global channel_bit_rate;
global Lsw;
global Lcw;
global Lbt;

% global L;
% global Arr_w;
% global Arr_b;
global k;
global scale;

global timeslot;
global Arr_b;
global Arr_w;
global Phi;

Pc = y(1);
alpha = y(2);
Pf = y(3);
Wc = y(4);
W0 = y(5);
% Arr_w = y(6);
% Arr_b = y(7);

x = (alpha + (1-alpha)*alpha/k)/1;


% mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
% mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
% mean_b = 3*Lbt - timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(alpha/x - (2*alpha/2*(alpha - 1))/x)*(x - 1) - (timeslot*(3*x - 3)*(Wi - 1)*(x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - (timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(3*x - 3)*(Wc - 1))/2;         
% mean_b = 3*Lbt - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - timeslot*(alpha/x - (2*alpha/k*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);

Srv_w = 1/mean_w;
Srv_b = 1/mean_b;
%     
%     Srv_w = 1/(mean_w)/N_w;
%     Srv_b = 1/(mean_b)/N_b;

Pw=1-min((Arr_w/100000)/Srv_w,1);
Pb=1-min((Arr_b/100000)/Srv_b,1);


b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw)*Pc/(1-Pc));
tau_w = Lsw*b001/(1-Pc);
% tau_w = b001/(1-Pc);

b000 = 1 / (scale*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ scale*Lbt+scale*(Losb-1)+scale/(1-Pb));%
tau_b = 3*Lbt*b000;


Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));

bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+max(Lcw,scale*Lbt)*Picbw);

Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));

BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+scale*Lbt*Bicb+scale*Lbt*Bisb+max(Lcw,scale*Lbt)*Bicbw);

eq1 = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b))) - Pc;
eq2 = (1-BMACii) - alpha;
eq3 = (1-bii) - Pf;
eq4 = tau_b/((1-tau_b)/Phi+tau_b) - tau_w;

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

%     disp(['  StdVar_w = ', num2str(sqrt(Var_Sw)) ]);%*N_w^2
    
    theta_w = (Arr_w/100000) - Srv_w;
%     delta_w = (Arr_w/100000)^3*(Arr_w/100000) + (Arr_w/100000)*Srv_w^2*Var_Sw; %Poisson
    delta_w = 0 + (Arr_w/100000)*Srv_w^2*Var_Sw; %Periodic *N_w^2
    
%     disp(['  theta_w = ', num2str((theta_w)) ]);
%     disp(['  delta_w = ', num2str(sqrt(delta_w)) ]);
    
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
%     disp(['  Q_w = ', num2str(Q_w) ]);
%     if(Q_w<0) 
%         Q_w=-1;
%     end
    
    delay_w = Q_w/(Arr_w/100000);
%     disp(['  delay_w = ', num2str(delay_w*10), ' us'  ]);

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

%     disp(['  StdVar_b = ', num2str(sqrt(Var_Sb)) ]);%*N_b^2
    
    
    theta_b = (Arr_b/100000) - Srv_b;
%     delta_b = (Arr_b/100000)^3*(Arr_b/100000) + (Arr_b/100000)*Srv_b^2*Var_Sb;%Poisson
    delta_b = 0 + (Arr_b/100000)*Srv_b^2*Var_Sb;%Periodic *N_b^2
    
%     disp(['  theta_b = ', num2str((theta_b)) ]);
%     disp(['  delta_b = ', num2str(sqrt(delta_b)) ]);
    
    
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
%     disp(['  Q_b = ', num2str(Q_b) ]);
%     if(Q_b<0) 
%         Q_b=-1;
%     end
    
    delay_b = Q_b/(Arr_b/100000);
%     disp(['  delay_b = ', num2str(delay_b*10), ' us' ]);




ceq = [
eq1;eq2;eq3;eq4%;eq5;eq6;eq7%
];
c = [-delay_b; -delay_w; -W0; -Wc];%



% ==================old===================
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
% 
% 
% x = alpha + (1-alpha)*alpha;
% b000 = 1 / (scale*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ scale*Lbt+scale*(Losb-1)+scale/(1-Pb));%
% eq1 = scale*Lbt*b000 - tau_b;
% 
% 
% Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
% Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
% Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
% Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
% Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));
% 
% BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+scale*Lbt*Bicb+scale*Lbt*Bisb+max(Lcw,scale*Lbt)*Bicbw);
% eq2 = (1-BMACii) - alpha;
% 
% 
% 
% b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw)*Pc/(1-Pc));
% eq4 = Lsw*b001/(1-Pc) - tau_w;
% 
% Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
% Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
% Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
% Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
% Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));
% 
% bii = 1/(1+Lsw*Pisw+Lcw*Picw+scale*Lbt*Picb+scale*Lbt*Pisb+max(Lcw,scale*Lbt)*Picbw);
% eq6 = (1-bii) - Pf;
% eq7 = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b))) - Pc;
% 
% 
% 
% timeslot = 1;
% 
% mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
% mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
%          
% 
% Srv_w = 1/(mean_w);
% Srv_b = 1/(mean_b);
% 
% eq5 = 1 - min((Arr_w/100000)/Srv_w,1) - Pw;
% eq3 = 1 - min((Arr_b/100000)/Srv_b,1) - Pb;
% 
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
% %     disp(['  StdVar_b = ', num2str(sqrt(Var_Sb)) ]);%*N_b^2
%     
%     
%     theta_b = (Arr_b/100000) - Srv_b;
% %     delta_b = (Arr_b/100000)^3*(Arr_b/100000) + (Arr_b/100000)*Srv_b^2*Var_Sb;%Poisson
%     delta_b = 0 + (Arr_b/100000)*Srv_b^2*Var_Sb;%Periodic *N_b^2
%     
% %     disp(['  theta_b = ', num2str((theta_b)) ]);
% %     disp(['  delta_b = ', num2str(sqrt(delta_b)) ]);
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
% %     disp(['  Q_b = ', num2str(Q_b) ]);
% %     if(Q_b<0) 
% %         Q_b=-1;
% %     end
%     
%     delay_b = Q_b/(Arr_b/100000);
% 
% 
% 
% 
% % Pw=0;
% % Pb=0;
% 
% %     Exp_Xb =  - (timeslot*x*(alpha/x - (2*alpha*(alpha - 1))/x))/(x - 1) + (timeslot*(3*x - 3)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(3*x - 3)*(Wc - 1))/(2*(x - 1)^2);
% % 
% %     Exp_Cb = alpha/x*1*timeslot+(1-alpha)*alpha/x*2*timeslot;
% % %     Var_Cb = alpha/x*(timeslot-Exp_Cb)^2+(1-alpha)*alpha/x*(2*timeslot-Exp_Cb)^2;
% % 
% %     ab = 3*timeslot*(Wi-1)/2;
% %     bb = Exp_Cb+3*timeslot*(Wc-1)/2;
% %     cb = Exp_Xb;
% % 
% %     Var_Sb = (3*timeslot^2*(Wi^2 - 1))/4 - (x*((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x - ...,
% %         (alpha*(alpha - 1)*(2*timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x))/(x - 1) - ...,
% %         (3*timeslot^2*x*(Wc^2 - 1))/(4*(x - 1)) + (ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1);
% % 
% % %     disp(['  StdVar_b = ', num2str(sqrt(Var_Sb)) ]);
% %     
% %     
% %     theta_b = Arr_b - Srv_b;
% %     delta_b = Arr_b^3*Arr_b + Arr_b*Srv_b^2*Var_Sb;
% %     
% % %     disp(['  theta_b = ', num2str((theta_b)) ]);
% % %     disp(['  delta_b = ', num2str(sqrt(delta_b)) ]);
% % 
% %     a = 2*theta_b/delta_b;
% %     f = @(z) -a.*z.*exp(a.*z);
% %     Q = integral(f,0,1000);
% % 
% %     delay = Q/Arr_b;
%     
%     
%     
% 
% 
% %=========method old===========
% 
% 
% 
% % tau_w = b001/(1-Pc);
% % tau_b = 3*b000;
% 
% 
%     ceq = [
%     eq1;eq2;eq3;eq4;eq5;eq6;eq7%
%     ];
%     c = [delay_b - 1000000];
%     
% 
% % a = x(1)*x(2)-x(1)-x(2);
% % b = -x(1)*x(2)-10;
% % 
% % 
% % % c = [ 
% % %     1.5+x(1)*x(2)-x(1)-x(2);
% % %     -x(1)*x(2)-10
% % %     ];
% % 
% % c = [
% %     1.5 + a;
% %     b
% %     ];
% % 
% % ceq = [];

