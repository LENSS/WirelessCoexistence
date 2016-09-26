
function f=coexist_math(xvect)
%suppose there are 2 wifi nodes and 6 bmac nodes.

m0=5;
W0=2^m0;
m=10;
b80211=0.826;
bbmac=0.174;
Lcw=10;
Lsw=8;
Low=10;
Ldifs=3;
Lsb=128*3;
Lob=210*3;
Wi=320;
Wc=80;
N=8;



alpha = xvect(1);
beta = xvect(2);
Px = xvect(3);

x = alpha + (1-alpha)*beta;
%b001=b80211 / (W0/2*((1-(2*Pc)^m))/(1-2*Pc)+(2*Pc)^m/2*(1-Pc)+Lcw*Pc/(1-Pc)+Lsw+Low +Ldifs*(1+Pc)/(1-Pc));
Backoff_terms = ((Px*W0*(1-(2*Px)^(m-1)))/(1-(2*Px))) + ( (Px + (W0*(2*Px)^m) + (2*((1-alpha)^Ldifs-Px))) / (2*(1-Px)) );
%DIFS_terms = ((((1 -(1-alpha)^Ldifs)*(1+Px))/(alpha*(1-Px))) );
DIFS_terms = (1+Px)/(1-Px) * (3-(3*alpha)+alpha^2);
Collision_terms = (Lcw*Px/(1-Px));
SuccessfulTx_terms = (Lsw*((1-alpha)^Ldifs-Px))/(1-Px);
OSDelay_terms = (Low*((1-alpha)^Ldifs-Px))/(1-Px);
b001 = b80211 / ( Backoff_terms + DIFS_terms + Collision_terms + SuccessfulTx_terms + OSDelay_terms);
b000=bbmac/(3*((Wi+1)/2+(x*(Wc+1))/(2*(1-x))+Lsb+Lob)); 
tau=b001*(1/(1-Px)) + 3*b000; 
f(1)=(1-(1-tau)^(N-1))*((Px*Lcw/(1-Px)+ Lsw*((1-alpha)^Ldifs-Px)/(1-Px))*b001 +3*Lsb*b000)-alpha; 
f(2)=(1-(1-tau)^N)/(2-(1-tau)^N)-beta;
f(3)=((1 - alpha)^(Ldifs) *(1-(1-tau)^(N-1)))-Px;


















