
function f=coexist_channel_math(xvect)
%suppose there are 2 wifi nodes and 6 bmac nodes.

% global SIFS;
% global PROP;
% global ACK;
% global ACKTO;
% global DIFS;
% 
global Wi;
global Wc;
% global Epl_b;
% global Epl_w;
global W0;
global m;
global N_w;
global N_b;
global Losw;
global Losb;
% global channel_bit_rate;

global Lsw;
global Lcw;
global Lbt;

global L;

Pc = xvect(1);
alpha = xvect(2);
Pf = xvect(3);

% Lbt = ceil(Epl_b*8*4/30)+PROP;
% Lsw = ceil(Epl_w*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACK;
% Lcw = ceil(Epl_w*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACKTO;


b001 = 1 /((W0)*(1-(2*Pc)^m)/(2*(1-2*Pc)*(1-Pf))+L*(W0*((2*Pc)^m)+1)/(2*(1-Pc)*(1-Pf))+(Lsw+Losw)+(Lcw+Losw)*Pc/(1-Pc));
tau_w = b001/(1-Pc);

x = alpha + (1-alpha)*alpha;
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 3*(1-alpha)/(1-x)+3*Lbt+3*Losb);%


%=========method old===========
tau_b = 3*b000;

%Ptrw = Lsw*(1-Pf)*b001+Lcw*(1-Pf)*Pc*b001/(1-Pc);%(Lsw*b001*(1-Pf) + Lcw*Pc*(1-Pf)*b001/(1-Pc))/(Lsw+Lcw);
%Ptrb = 3*Lbt*b000;
if N_b == 0
    Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1));
    Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)));

    bii = 1/(1+Lsw*Pisw+Lcw*Picw);
    f(1) = (1 - ((1-tau_w)^(N_w-1))) - Pc;
    f(3) = (1-bii) - Pf;
    
elseif N_w == 0
    Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1));
    Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)));

    BMACii = 1/(1+3*Lbt*Bicb+3*Lbt*Bisb);
    f(1) = (1 - ((1-tau_b)^(N_b))) - Pc;
    f(2) = (1-BMACii) - alpha;
    
else
    Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*((1-tau_b)^(N_b));
    Pisb = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^(N_w-1));
    Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*((1-tau_b)^(N_b));
    Picb = (1 - ((1-tau_b)^N_b) - N_b*tau_b*((1-tau_b)^(N_b-1)))*((1-tau_w)^(N_w-1));
    Picbw = (1 - ((1-tau_b)^(N_b)))*(1-((1-tau_w)^(N_w-1)));

    bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);

    Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b-1));
    Bisb = (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1))*((1-tau_w)^(N_w));
    Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*((1-tau_b)^(N_b-1));
    Bicb = (1 - ((1-tau_b)^(N_b-1)) - (N_b-1)*tau_b*((1-tau_b)^(N_b-1-1)))*((1-tau_w)^(N_w));
    Bicbw = (1 - ((1-tau_b)^(N_b-1)))*(1-((1-tau_w)^(N_w)));

    BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);

    f(1) = (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b))) - Pc;
    f(2) = (1-BMACii) - alpha;
    f(3) = (1-bii) - Pf;
end
%=========method new===========

% % % P_I = 1 - alpha;
% % % phi = 3*b000/(1-x);
% % % 
% % % Pisw = (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b)));
% % % Pisb = 3*Lbt*N_b*phi*P_I^2*(1-phi)^(N_b-1)*((1-tau_w)^(N_w-1));
% % % Picw = (1 - ((1-tau_w)^(N_w-1)) - (N_w-1)*tau_w*((1-tau_w)^(N_w-1-1)))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b)));
% % % Picb = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b))) - 3*Lbt*N_b*phi*P_I^2*(1-phi)^(N_b-1))*((1-tau_w)^(N_w-1));
% % % Picbw = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b))))*(1-((1-tau_w)^(N_w-1)));
% % % 
% % % bii = 1/(1+Lsw*Pisw+Lcw*Picw+3*Lbt*Picb+3*Lbt*Pisb+3*Lbt*Picbw);
% % % 
% % % Bisw = N_w*tau_w*((1-tau_w)^(N_w-1))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b-1)));
% % % Bisb = 3*Lbt*(N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1)*((1-tau_w)^(N_w));
% % % Bicw = (1 - ((1-tau_w)^N_w) - N_w*tau_w*((1-tau_w)^(N_w-1)))*(1-3*Lbt*P_I^2*(1-(1-phi)^(N_b-1)));
% % % Bicb = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b-1))) - 3*Lbt*(N_b-1)*phi*P_I^2*(1-phi)^(N_b-1-1))*((1-tau_w)^(N_w));
% % % Bicbw = (1 - (1-3*Lbt*P_I^2*(1-(1-phi)^(N_b-1))))*(1-((1-tau_w)^(N_w)));
% % % 
% % % BMACii = 1/(1+Lsw*Bisw+Lcw*Bicw+3*Lbt*Bicb+3*Lbt*Bisb+3*Lbt*Bicbw);
% % % 
% % % f(1) = (1 - ((1-tau_w)^(N_w-1))*((1-3*Lbt*P_I^2*(1-(1-phi)^(N_b))))) - Pc;

%=================================

% f(1) =  (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b))) - Pc;%(1 - ((1-Ptrw)^(N_w-1)))+(1-((1-Ptrb)^(N_b)))-
% f(2) =  (1 - ((1-Ptrw)^(N_w))*((1-Ptrb)^(N_b-1))) - alpha;% (1 - ((1-Ptrw)^N_w))+(1-((1-Ptrb)^(N_b)))-
% f(3) =  (1 - ((1-Ptrw)^(N_w-1))*((1-Ptrb)^(N_b)))  - Pf;%(1 - ((1-Ptrw)^(N_w)))+(1-((1-Ptrb)^(N_b)))-
% % % f(1) = (Lcw*Picw+3*Lbt*Picbw)*bii - Pc;
end
