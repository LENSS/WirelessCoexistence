clear all;

m0=5;
W0=2^m0;
m=10;
n=15;
Wi=320;
Wc=80;
%Lsw=30;
%Lcw=36;
Losw=20;
%Lsw=1;
%Lcw=1;
%Losw=0;
%Lbt=128;
Losb=210;
% N_w=1;
N_w=2;
Trials = 4;

averThr_w(1:Trials)=0;
aggrThr_w(1:Trials)=0;
averThr_b(1:Trials)=0;
aggrThr_b(1:Trials)=0;
alphas(1:Trials)=0;
tau_ws(1:Trials)=0;
tau_bs(1:Trials)=0;

for N_b=1:Trials

xguess=[0 0 0 N_w N_b]';
xvect = fsolve('coexist_math_simplest', xguess);
alpha = xvect(2);
Pc = xvect(1);
Pf = xvect(3);
disp('======================================================================');
disp('======================================================================');
disp(['  N_w = ', num2str(N_w) ])
disp(['  N_b = ', num2str(N_b) ])

disp('The roots from the default "fsolve" are: ')
disp(['  alpha = ', num2str(alpha) ])
disp(['  Pc = ', num2str(Pc) ])
disp(['  Pf = ', num2str(Pf) ])

Lbt = ceil(Epl_b*8*4/30);
Lsw = ceil(Epl_w*8/54/10);
Lcw = ceil(Epl_w*8/54/10);

x = alpha + (1-alpha)*alpha;
b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1 / ((W0/2)*(1-(2*Pc)^m)/(1-2*Pc)+(W0*((2*Pc)^m)+1)/(2*(1-Pc))+(Lsw+Losw)*(1-Pf) +(Lcw*Pc)*(1-Pf)/(1-Pc));
%b001 = 1/ ((W0*(1-((2*Pc)^m+1))/(2*(1-2*Pc)) ) + ( (2^(m+1))*W0*(1-(Pc^(n-m)))*(Pc^(m+1))/(1-Pc) ) + ((1-(Pc^(n+1)))/(2*(1-Pc))) + ((1-(Pc^(n+1)))*(1-Pf)*(Lsw+(Lcw*Pc))/(1-Pc)) );
%tau_w = b001*(1-Pf)*(1-(Pc^(n+1)))/(1-Pc);
tau_w = b001*(1-Pf)/(1-Pc);
b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+ 4*(1-alpha)/(1-x)+3*Lbt+3*Losb);
%b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+3*Lbt+2+3*Losb);
%b000 = 1 / (3*((Wi+1)/2+(Wc+1)*x/(2*(1-x))+(1-alpha)/(1-x))+3*Lbt+2+3*Losb);
tau_b = b000;%/(1-x);

disp(['  tau_w = ', num2str(tau_w) ])
disp(['  tau_b = ', num2str(tau_b) ])

Ptrw = Lsw*b001*(1-Pf)+Lcw*Pc*b001*(1-Pf)/(1-Pc);%(Lsw*b001*(1-Pf) + Lcw*Pc*(1-Pf)*b001/(1-Pc))/(Lsw+Lcw);
Ptrb = 3*Lbt*b000; 
disp(['  Ptrw = ', num2str(Ptrw) ])
disp(['  Ptrb = ', num2str(Ptrb) ])

alphas(N_b) = alpha;
tau_ws(N_b) = tau_w;
tau_bs(N_b) = tau_b;

Ptr_w = 1-((1-tau_w)^N_w);%*(1-tau_b)^N_b;
Psucc_w = N_w*tau_w*((1-tau_w)^(N_w-1))*((1-tau_b)^N_b)/Ptr_w;
Pcolli_w = 1-Psucc_w;
Posd_w = Psucc_w+Pcolli_w;
Pbo_w = (1-tau_w)^N_w;

disp(['  Ptr_w = ', num2str(Ptr_w) ])
disp(['  Psucc_w = ', num2str(Psucc_w) ])
disp(['  Pcolli_w = ', num2str(Pcolli_w) ])
disp(['  Posd_w = ', num2str(Posd_w) ])
disp(['  Pbo_w = ', num2str(Pbo_w) ])

Ts_w = 280;
SIFS = 10;
PROP = 1;
ACK = 20;
%ACKTO = 40;
Ts_w_max = 1120;%??? why 1280.
DIFS = 30;
%Tsucc_w = Ts_w+SIFS+PROP+ACK+DIFS+PROP;
%Tcolli_w = Ts_w_max+DIFS+PROP;%??? ACKTIMEOUT?
%Tbo_w = 10;
Tsucc_w = 350;
Tcolli_w = 350;
Tbo_w = 10;
Tosd_w = 380;
Epl_w = 2346;

Throughput_w = Ptr_w*Psucc_w*Epl_w/(Ptr_w*Psucc_w*Tsucc_w+Ptr_w*Pcolli_w*Tcolli_w+Pbo_w*Tbo_w+Ptr_w*Posd_w*Tosd_w);
Max_tput_w = Ptr_w*Epl_w*8/(Ptr_w*Tsucc_w+Pbo_w*Tbo_w+Ptr_w*Posd_w*Tosd_w);
averThr_w(N_b) = Throughput_w*8/N_w/54;%/Max_tput_w;
aggrThr_w(N_b) = Throughput_w*8/54;%/Max_tput_w;
disp(['  Max_throughput_w = ', num2str(Max_tput_w) , ' Mbps'])
disp(['  Throughput_w = ', num2str(Throughput_w*8) , ' Mbps'])
disp(['  Throughput_w/N = ', num2str(Throughput_w*8/N_w) , ' Mbps'])

Ptr_b = 1-(1-tau_b)^N_b;%*(1-tau_w)^N_w;
Psucc_b = N_b*tau_b*((1-tau_b)^(N_b-1))*((1-tau_w)^N_w)/Ptr_b;
Pcolli_b = 1-Psucc_b; %???
Posd_b = Psucc_b+Pcolli_b; %???
Pbo_b = (1-tau_b)^N_b;

disp(['  Ptr_b = ', num2str(Ptr_b) ])
disp(['  Psucc_b = ', num2str(Psucc_b) ])
disp(['  Pcolli_b = ', num2str(Pcolli_b) ])
disp(['  Posd_b = ', num2str(Posd_b) ])
disp(['  Pbo_b = ', num2str(Pbo_b) ])

Ts_b = 136*30;
Tsucc_b = Ts_b;
Tcolli_b = Ts_b; %???
Tbo_b = 30;
Tosd_b = 220*30;
Epl_b = 128;

Throughput_b = Ptr_b*Psucc_b*Epl_b/(Ptr_b*Psucc_b*Tsucc_b+Ptr_b*Pcolli_b*Tcolli_b+Pbo_b*Tbo_b+Ptr_b*Posd_b*Tosd_b);
Max_tput_b = Ptr_b*Epl_b*8000/(Ptr_b*Tsucc_b+Pbo_b*Tbo_b+Ptr_b*Posd_b*Tosd_b);
averThr_b(N_b) = Throughput_b*8*1000/N_b/250;%/Max_tput_b;
aggrThr_b(N_b) = Throughput_b*8*1000/250;%/Max_tput_b;
disp(['  Max_throughput_b = ', num2str(Max_tput_b) , ' Kbps'])
disp(['  Throughput_b = ', num2str(Throughput_b*8*1000) , ' Kbps'])
disp(['  Throughput_b/N = ', num2str(Throughput_b*8*1000/N_b) , ' Kbps'])
disp('======================================================================');
end

plot(aggrThr_w,'DisplayName','aggrThr_w','YDataSource','aggrThr_w');hold all;
%plot(averThr_w,'DisplayName','averThr_w','YDataSource','averThr_w');hold all;
plot(aggrThr_b,'DisplayName','aggrThr_b','YDataSource','aggrThr_b');hold all;
%plot(averThr_b,'DisplayName','averThr_b','YDataSource','averThr_b');
plot(alphas,'DisplayName','alpha','YDataSource','alphas');hold all;
plot(tau_ws,'DisplayName','tau_w','YDataSource','tau_ws');hold all;
plot(tau_bs,'DisplayName','tau_b','YDataSource','tau_bs');hold all;
hold off;grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,50,0,1]);figure(gcf);