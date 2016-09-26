
function f=coexist_math_retry_limit(xvect)
%suppose there are 2 wifi nodes and 6 bmac nodes.

m0=5;
W0=2^m0;
m=10;
n=15;
% b80211=0.826;
% bbmac=0.174;
% Lcw=10;
% Lsw=8;
% Low=10;
% Ldifs=3;
% Lsb=128*3;
% Lob=210*3;
Wi=320;
Wc=80;
N_w=1;
N_b=4;

alpha = xvect(1);
Pc = xvect(2);
Pf = xvect(3);

x = alpha + (1-alpha)*alpha;
b001 = 1 / ( (W0*(1-((2*Pc)^m+1))/(2*(1 - 2*Pc)) ) + ( (2^(m+1))*W0*(1-(Pc^(n-m)))*(Pc^(m+1))/(1-Pc) ) + ((1-(Pc^(n+1)))/(2*(1-Pc))) + ((1-(Pc^(n+1)))*(1-Pf)*(1+Pc)/(1-Pc) ) );
tau_w = b001/(1-Pc);
b000 = 1 / ((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x)+1);
tau_b = b000/(1-x);

f(1) = 1 - ((1-tau_w)^(N_w-1)*(1-tau_b)^N_b) - Pc;
f(2) = 1 - (((1-tau_w)^N_w)*(1-tau_b)^(N_b-1)) - alpha;
f(3) = 1 - ((1-tau_w)^(N_w-1)*(1-tau_b)^N_b) - Pf;












