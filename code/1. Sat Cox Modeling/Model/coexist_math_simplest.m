
function f=coexist_math_simplest(xvect)
%suppose there are 2 wifi nodes and 6 bmac nodes.

%m0=4;
%W0=2^m0;
%m=6;
n=15;
% b80211=0.826;
% bbmac=0.174;
% Lcw=10;
% Lsw=8;
% Low=10;
% Ldifs=3;
% Lsb=128*3;
% Lob=210*3;
%Wi=320;
%Wc=80;
% N_w=4;
% N_b=0;
%Lsw=30;
%Lcw=36;
% Losw=60;
%Lsw=1;
%Lcw=1;
%Losw=0;
%Lbt=128;
% Losb=220;
channel_bit_rate=54;
Pc = xvect(1);
alpha = xvect(2);
Pf = xvect(3);
Wi = xvect(4);
Wc = xvect(5);
Epl_b = xvect(6);
Epl_w = xvect(7);
Losw = xvect(12);
Losb = xvect(13);

m0=xvect(8);
W0=2^m0;
m=xvect(9);

N_w = xvect(10);
N_b = xvect(11);

Lbt = 1;%ceil(Epl_b*8*4/30);
Lsw = 1;%ceil(Epl_w*8/channel_bit_rate/10);
Lcw = 1;%ceil(Epl_w*8/channel_bit_rate/10);

x = alpha + (1-alpha)*alpha;
b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1/ ((W0*(1-((2*Pc)^m+1))/(2*(1-2*Pc)) ) + ( (2^(m+1))*W0*(1-(Pc^(n-m)))*(Pc^(m+1))/(1-Pc) ) + ((1-(Pc^(n+1)))/(2*(1-Pc))) + ((1-(Pc^(n+1)))*(1-Pf)*(Lsw+(Lcw*Pc))/(1-Pc)) );
%tau_w = b001*(1-Pf)*(1-(Pc^(n+1)))/(1-Pc);
tau_w = b001*(1-Pf)/(1-Pc);
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 4*(1-alpha)/(1-x)+3*Lbt+3*Losb);%
%b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+3*Lbt+2+3*Losb);
tau_b = b000;%/(1-x);
%phi_b = b000/(1-x);
Ptrw = Lsw*(1-Pf)*b001+Lcw*(1-Pf)*Pc*b001/(1-Pc);%(Lsw*b001*(1-Pf) + Lcw*Pc*(1-Pf)*b001/(1-Pc))/(Lsw+Lcw);
Ptrb = 3*Lbt*b000;
f(1) =  (1 - ((1-Ptrw)^(N_w-1))*((1-Ptrb)^(N_b))) - Pc;%(1 - ((1-Ptrw)^(N_w-1)))+(1-((1-Ptrb)^(N_b)))-
f(2) =  (1 - ((1-Ptrw)^(N_w))*((1-Ptrb)^(N_b))) - alpha;% (1 - ((1-Ptrw)^N_w))+(1-((1-Ptrb)^(N_b)))-
f(3) =  (1 - ((1-Ptrw)^(N_w))*((1-Ptrb)^(N_b)))  - Pf;%(1 - ((1-Ptrw)^(N_w)))+(1-((1-Ptrb)^(N_b)))-












