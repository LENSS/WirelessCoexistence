function f = delay_objfun(y)

global Wi;

global m;
global N_w;
global N_b;
global Losw;
global Losb;
% global channel_bit_rate;
global Lsw;
global Lcw;
global Lbt;
global timeslot;
global k;
global scale;
global Arr_b;

Pc = y(1);
alpha = y(2);
Pf = y(3);
Wc = y(4);
W0 = y(5);
Arr_w = y(6);
% Arr_b = y(7);


x = (alpha + (1-alpha)*alpha/k)/1;

mean_w = Lsw - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)));
mean_b = scale*Lbt - (timeslot*x*(alpha/x - (2*alpha/k*(alpha - 1))/x))/(x - 1) + (timeslot*(scale*x - scale)*(Wi - 1))/(2*(x - 1)) - (timeslot*x*(scale*x - scale)*(Wc - 1))/(2*(x - 1)^2);
% mean_b = 3*Lbt - timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(alpha/x - (2*alpha/2*(alpha - 1))/x)*(x - 1) - (timeslot*(3*x - 3)*(Wi - 1)*(x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - (timeslot*(x/(x - 1)^2 + (x^11*(10*x - 11))/(x - 1)^2)*(3*x - 3)*(Wc - 1))/2;         
% mean_b = 3*Lbt - (timeslot*(3*x - 3)*(Wc - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x))/2 - (timeslot*(3*x - 3)*(Wi - 1)*(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1))/2 - timeslot*(alpha/x - (2*alpha/k*(alpha - 1))/x)*(x - 1)*(6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x);

Srv_w = 1/(mean_w);
Srv_b = 1/(mean_b);
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


Pi = (1-tau_w)^(N_w)*(1-tau_b)^(N_b);
Psw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
Psb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w));
Pcw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b));
Pcb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w));
Pcbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w)));


Psucc_w = Psw/(Psw*Lsw+Psb*scale*Lbt+Pcw*Lcw+Pcb*scale*Lbt+Pcbw*max(Lcw,scale*Lbt)+Pi);


f = -Psucc_w;


% =============old=============
% tau_w = y(1); %tau_w tau_b alpha Pb Wc
% tau_b = y(2);
% 
% 
% % 
% % 
% % b001 = 1/((W(1))*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W(1)*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw-1)+1/(1-Pw)+(Lcw+Losw)*Pc/(1-Pc));
% % b000 = 1/(3*((Wi+1)/2+(W(2)+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*(Losb-1)+3/(1-Pb));%
% % 
% % 
% % 
% % 
% f = -N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b));
% 


% t = (4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
% 
% f=exp(x(1))*t;
