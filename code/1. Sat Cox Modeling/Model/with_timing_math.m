
function f=with_timing_math(xvect)
%suppose there are 2 wifi nodes and 6 bmac nodes.

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

SIFS = 1;
PROP = 1;
ACK = 2;
ACKTO = 50;
DIFS = 6;

Lbt = ceil(Epl_b*8*4/30);
Lsw = ceil(Epl_w*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACK;
Lcw = ceil(Epl_w*8/channel_bit_rate/10)+SIFS+PROP+DIFS+ACKTO;

x = alpha + (1-alpha)*alpha;

b001 = 1 /((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));

tau_w = b001*(1-Pf)/(1-Pc);
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 7*(1-alpha)/(1-x)+3*Lbt+3*Losb);%

tau_b = b000;

Ptrw = Lsw*(1-Pf)*b001+Lcw*(1-Pf)*Pc*b001/(1-Pc);%(Lsw*b001*(1-Pf) + Lcw*Pc*(1-Pf)*b001/(1-Pc))/(Lsw+Lcw);
Ptrb = 3*Lbt*b000;

f(1) =  (1 - ((1-tau_w)^(N_w-1))*((1-tau_b)^(N_b))) - Pc;%(1 - ((1-Ptrw)^(N_w-1)))+(1-((1-Ptrb)^(N_b)))-
f(2) =  (1 - ((1-Ptrw)^(N_w))*((1-Ptrb)^(N_b-1))) - alpha;% (1 - ((1-Ptrw)^N_w))+(1-((1-Ptrb)^(N_b)))-
f(3) =  (1 - ((1-Ptrw)^(N_w-1))*((1-Ptrb)^(N_b)))  - Pf;%(1 - ((1-Ptrw)^(N_w)))+(1-((1-Ptrb)^(N_b)))-
end












