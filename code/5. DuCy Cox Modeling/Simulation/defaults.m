
%defaults:
L_W = 15+2;% WiFi time slots (15 == 1000 bytes)
Ackto_W = 4;% WiFi time slots

cw_W = 16;
CW_W(1:2) = [cw_W 16];%minimal CW size

N_W = 40; % N_W/2 senders
T_W = 100;% ms
Rho_W = 0.1;
TR_W = 5;

L_B = 34; %10;%100;%BMAC time slots (100 == 94 bytes)
Ackto_B = 10;%256; % BMAC time slots (= 300us)

cw_B = 70;
CW_B(1:2) = [cw_B 70];

N_B = 40; % N_W/2 senders
T_B = 400;% ms
Rho_B = 0.05;
TR_B = 1;


round = 5000; %set round to simulate
sim_time = round*T_W;%10000;%ms