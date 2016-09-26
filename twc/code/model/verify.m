%function F = myfun(x)
%F = [2*x(1) - x(2) - exp(-x(1));
%      -x(1) + 2*x(2) - exp(-x(2))];

% 
% function f=myfun(xvect)
% x = xvect(1);
% y = xvect(2);
% z = xvect(3);
% xyz = x + y + z
% f(1)=sin(x) + y^2 + log(z) - 7 
% f(2)=3*x + 2^y - z^3 + 1     
% f(3)=xyz - 5


% function F = myfun(x,k)
% H=0.32;
% Pc0=0.23;W=0.18;
% F=[Pc0+H*(1+1.5*(x(1)/W-1)-0.5*(x(1)/W-1)^3)-x(2);
%  x(1)-k*sqrt(x(2))]

function f=verify(xvect)
m = 4;
%mb = 8;
m0 = 3;
n = 3;
L = 7;
N = 10;
Ls = L + 6;
Lc = L + 4;
Lack = 2;
W0 = 2^m0;
alpha = xvect(1);
beta = xvect(2);
Pc = xvect(3);
x = alpha + (1-alpha)*beta;
y = Pc*(1-x^(m+1));
temp1 = ((1-(2*x)^(m+1))/(1-2*x)*W0+(1-x^(m+1))/(1-x))*(1-y^(n+1))/(1-y);
temp2 = (1-alpha)*((1-x^(m+1))/(1-x))*((1-y^(n+1))/(1-y));
temp3 = ((1-y^(n+1))/(1-y))*(1-x^(m+1))*(Ls*(1-Pc)+Lc*Pc);
b000 = 1/(temp1/2+temp2+temp3);
tau = b000*((1-x^(m+1))/(1-x))*((1-y^(n+1))/(1-y));
f(1) = 1-(1-tau)^(N-1)-Pc;
alpha1 = L*(1-(1-tau)^(N-1))*(1-alpha)*(1-beta);
alpha2 = Lack*(N*tau*((1-tau)^(N-1))/(1-(1-tau)^N))*(1-(1-tau)^(N-1))*(1-alpha)*(1-beta);
f(2) = alpha1+alpha2-alpha;
f(3) = beta-(Pc+N*tau*((1-tau)^(N-1)))/(2-(1-tau)^N+(N*tau*((1-tau)^(N-1))));














