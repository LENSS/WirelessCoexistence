
clear all;
clear all;
clear all;

% globals;
% global Service;

ii = 0;
jj = 0;

% no_B = 20; % 10:5:30% 5:5:40 % %number of senders
% to_B = 300; %timeout, us
% cy_B = 400*1e3; %cycle length, us
% ra_B = 0.1; %duty cycle ratio
% cw_B = 70; %CW size
% Ds_B = 1020;%3000;%us

% no_W = 20; %number of senders
% to_W = 10; %timeout, us
% cy_W = 100*1e3; %cycle length, us
% ra_W = 0.1; %duty cycle ratio
% cw_W = 16; %minimal CW size

% for N = 20:10:20 %# senders
% 
% Pe(1:N) = 0;
% P_free(1:N) = 0;
% Service(1:N) = 0;
% S(1:N) = 0;


% t = 0;
% L_bmac = 100;%bytes
L_B = 1020;%3000;%us
Ackto_B = 300;%256; %us
% timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot
L_W = 150; %Packet size in DATA us
Ackto_W = 10;
sim_time = 2000000;%ms

% mac_queue(1:N) = 0; %# packets
% QUEUE_SIZE = 100;

global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
global NN_B NN_W;
global TR_B TR_W;
global q0b q0w;

% global actualusage_d;
% global energyusage_d LL_Wd;
global TUNE_FACTOR
TUNE_FACTOR = 800;%0;%
% global PERCENTAGE;
% PERCENTAGE = 0.5;

DRAW_ONLY = 1;

for T_W = 100000  %<<<<<<<<<<<<==================== cycle of PSM (us)
    if DRAW_ONLY
        break;
    end
    for Rho_W = 0.1 %[0.1 0.2] %duty cycle ratio
        for cw_W = 16
            CW_W(1:2) = [cw_W 16];
            for NN_W = 20 %5:5:30 %<<<<<<<<<<<<==================== % # senders
                for TR_W = 1:2:11 %0.001%7%
                    for T_B = 400000 % 100000:100000:400000 %<<<<<<<<<<<<==================== cycle of LPL (us) 
                        for Rho_B = 0.05 %[0.05 0.1] %duty cycle ratio
                            for cw_B = 70 %[70 310] %
                                CW_B(1:2) = [cw_B 70];

%                                 jj = jj + 1;
%                                 ii = 0;

                                for TR_B = [0.5 1 3] %0.5%[2.1 2.2 2.3 2.4 2.5] %[0.1 0.2] %
                                    for NN_B = 20%5:5:30 %<<<<<<<<<<<<==================== % # senders

%                                     Pe(1:N_B) = 0;
%                                     P_free(1:N_B) = 0;
%                                     Service(1:N_B) = 0;
%                                     Service2(1:N_B) = 0;
                                        fprintf('=============================================\n');
                                        fprintf('# LPL Nodes: %d (%d senders)\n', NN_B*2, NN_B);
                                        fprintf('CW: %d time slots\n', CW_B(1));
                                        fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
                                        fprintf('Traffic: %g p/s\n', TR_B);
                                        fprintf('# PSM Nodes: %d (%d senders)\n', NN_W*2, NN_W);
                                        fprintf('CWmin: %d time slots\n', CW_W(1));
                                        fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
                                        fprintf('Traffic: %g p/s\n', TR_W);
                                        fprintf('=============================================\n');

                                        outname=sprintf('result-(%g)-(%d-%g-%d)-(%d-%g-%d)-(%d).mat', ...
                                                                   TR_B, T_B, Rho_B, CW_B(1), ...
                                                                    T_W, Rho_W, CW_W(1), TUNE_FACTOR);
                                                                
                                        if 0%exist(outname, 'file')%
                                          % File exists.  Do stuff....
                                          load(outname);
                                        else
                                          % File does not exist.
                                            Sr_B(1:NN_W+1,1:NN_B+1) = 0;
                                            Pe_B(1:NN_W+1,1:NN_B+1) = 0;
                                            Sr_W(1:NN_W+1,1:NN_B+1) = 0;
                                            En_W(1:NN_W+1,1:NN_B+1) = 0;

                                        %     for traffic_rate_bmac = [0.2 0.5 1 2] % X = X packets/s
                                            for n_B = NN_B % 0:NN_B % 
                                                for n_W = NN_W % 0:NN_W % 
                                                    if  n_W~=0 || n_B~=0
                                                        out = model(sim_time, n_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                                                                n_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1));
                                                    end
                                                end
        %                                         P_free(n) = out(1);
        %                                         Pe(n) = out(2);
        %                                         Service2(n) = out(3);
        %                                         Service(n) = out(4);
                                        %         fprintf('beta %i: %g\n', n, nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n));
                                        %         SS = SS + nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n)*S(n);
                                            end
%                                             Pe_B(NN_W+1,NN_B+1)

                                            return;
                                            
                                            save(outname, 'Sr_B', 'Pe_B', 'Sr_W', 'En_W', 'P_Wa', 'P_Wd');

                                        
                                        end
                                            
                                        
%                                         close all;
%                                         end
                                        continue;
                                        
                                        Q0_B = 0.5;
                                        Q0_W = 0.5;
                                        
%                                         for i=0:NN_B
%                                             beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
%                                         end
% 
%                                         for j=0:NN_W
%                                             beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
%                                         end                                        
%                                         
%                                         temp = beta_W * Sr_B;
%                                         temp1 = beta_B * Sr_W';
%                                         
%                                         SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
%                                         SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
%                                         
%                                         continue;

                                        xguess=[Q0_B Q0_W]';
                                        options = optimoptions('fsolve','Display','off'); %'iter'

                                        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);
                                        
                                        if exitflag<=0
%                                            error('adfad');
%                                             ii = ii + 1;
%                                             result(ii,1) = 0;
%                                             result(ii,2) = 0;
%                                             result(ii,3) = 0;
%                                             result(ii,4) = 0;
%                                             zzww{ii,1} = sprintf('N_W:%d,TR_W:%d,N_B:%d,TR_B:%g', NN_W, TR_W, NN_B, TR_B);
%                                            continue;

                                            Q0_B = q0b;
                                            Q0_W = q0w;

                                        else

                                            Q0_B = xvect(1);
                                            Q0_W = xvect(2);
                                        end
                                        
                                        fprintf('Expected Q0_B: %g\n', Q0_B);
                                        fprintf('Expected Q0_W: %g\n', Q0_W);
%                                         NN_B
%                                         NN_W
                                        for i=0:NN_B
                                            beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
                                        end

                                        for j=0:NN_W
                                            beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
                                        end


                                        temp = beta_W(1:NN_W+1) * Pe_B(1:NN_W+1,1:NN_B+1);
%                                         size(beta_B)
%                                         size(Sr_W(1:NN_W+1,1:NN_B+1))
                                        temp1 = beta_B(1:NN_B+1) * (1e6./Sr_W(1:NN_W+1,1:NN_B+1))';
                                        temp2 = beta_B(1:NN_B+1) * (En_W(1:NN_W+1,1:NN_B+1))';
%                                         temp1 = beta_B * Sr_W(1:NN_W+1,1:NN_B+1)';

                                        SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
                                        SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
                                        en_W = beta_W(2:NN_W+1) * temp2(2:NN_W+1)';

%                                         beta_B(1)
%                                         1-SS_B
                                        SS_B = SS_B/(1-beta_B(1));
                                        SS_W = SS_W/(1-beta_W(1));
                                        en_W = en_W/(1-beta_W(1));
                                        
                                        TH_B(NN_W,NN_B) = (1-SS_B);
                                        
                                        fprintf('Average Pe_B: %g p/s\n', 1-TH_B(NN_W,NN_B));
                                        
                                        fprintf('Average throughput of LPL: %g p/s\n', TH_B(NN_W,NN_B)*TR_B);
                                        fprintf('Average throughput of LPL (percent): %g \n', TH_B(NN_W,NN_B));
                                       
                                        EN_B(NN_W,NN_B) = 1-Q0_B;
                                        
                                        fprintf('Average energy consumption of LPL: %g\n', EN_B(NN_W,NN_B));
                                        
                                        TH_W(NN_W,NN_B) = SS_W; %classical way to compute service rate (i.e. throughput)
%                                         TR_W/(1 - Q0_W) % the reciprocal of average service time
                                        
                                        fprintf('Average throughput of PSM: %g p/s\n', TH_W(NN_W,NN_B));
                                        
                                        
                                        EN_W(NN_W,NN_B) = en_W;
                                        fprintf('Average energy consumption of PSM: %g\n', EN_W(NN_W,NN_B));
                                        ii = ii + 1;
                                        result(ii,1) = TH_B(NN_W,NN_B);
                                        result(ii,2) = TH_W(NN_W,NN_B);
                                        result(ii,3) = EN_B(NN_W,NN_B);
                                        result(ii,4) = EN_W(NN_W,NN_B);
                                        zzww{ii,1} = sprintf('N_W:%d,TR_W:%d,N_B:%d,TR_B:%g', NN_W, TR_W, NN_B, TR_B);
                                        
%                                         En_B(NN_W,NN_B) = 1-Q0_B;
%                                         
%                                         fprintf('Average energy consumption of PSM: %g\n', En_B(NN_W,NN_B));



%                                         P0 = 0.01;
% 
%                                         xguess=[P0]';
%                                         options = optimoptions('fsolve','Display','off'); %'iter'
% 
%                                         xvect = fsolve('servtime_solver', xguess, options);
% 
%                                         P0 = xvect(1);
%                                         fprintf('Expected P0: %g\n', P0);
% 
%                                         ii = ii + 1;
% 
%                                         SS = 0;
%                                         for i=1:N_B
%                                           SS = SS + nchoosek(N_B,i)*(1-P0)^i*P0^(N_B-i)*Pe(i);
%                                         end
%                                         thr(ii,jj) = 1-SS;
%                                         fprintf('Expected S: %g\n', thr(ii,jj));
% 
%                     %                     SSS = 0;
%                     %                     for i=0:N-1
%                     %                       SSS = SSS + nchoosek(N-1,i)*(1-P0)^i*P0^(N-1-i)*Service2(i+1);
%                     %                     end
%                     % %                     ene(ii,jj) = (1-P0)*SSS/T;
%                     %                     ene(ii,jj) = SSS/T;
% 
%                     %                     fprintf('Expected Servetime: %g\n', ((1-SS)*T/2 + SS*T));
%                     %                     
%                     %                     ene(ii,jj) = (1-P0)*((1-SS)*T/2 + SS*T)/T;
%                     % %                     ene(ii,jj) = (1-P0)*2*((1-SS)*T/2 + SS*T)/T;
% 
%                                         ene(ii,jj) = (1-P0);
%                                         fprintf('Expected energy: %g\n', ene(ii,jj));
% 
%                     %                     E_Service =(1-P0)/traffic_rate_bmac;
%                     %                     fprintf('E_Service: %g\n', E_Service*1e6);
%                     %                     ene(ii,jj) = (1-P0)*E_Service/T*1e6;
%                     %                     fprintf('Expected energy: %g\n', ene(ii,jj));
% 
% 
% 



                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if DRAW_ONLY ~= 1
   return; 
end
% return;

% 
% clear all;
% clear all;
% clear all;
% 
% global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
% global NN_B NN_W;
% global TR_B TR_W;
% global T_B T_W;
% global q0b q0w;
% global CASE;
% global TUNE_FACTOR TUNE_FACTOR2;
% 
% %defaults:
% L_W = 150; %Packet size in DATA us
% Ackto_W = 10;
% 
% cw_W = 16;
% CW_W(1:2) = [cw_W 16];%minimal CW size
% 
% N_W = 40; % N_W/2 senders
% T_W = 100*1e3;% ms
% Rho_W = 0.1;
% TR_W = 5;
% 
% L_B = 1020;%3000;%us
% Ackto_B = 300;%256; %us
% 
% cw_B = 70;
% CW_B(1:2) = [cw_B 70];
% 
% N_B = 40; % N_W/2 senders
% T_B = 400*1e3;% ms
% Rho_B = 0.05;
% TR_B = 1;
% 
% TUNE_FACTOR = 800;%0;%
% TUNE_FACTOR2 = 250;%0;%
% 
% %Figure type 1 (vary # nodes)
% % x-axis: # PSM (5 10 15 20)
% % each line: # LPL (10 20 30)
% % y-axis: Throughput of PSM
% % y-axis: Energy consumption of PSM
% % y-axis: Throughput of LPL
% % y-axis: Energy consumption of LPL
% ii = 1;
% jj = 1;
% 
% result_ThrB{1,jj} = 'Thr_B';
% result_ThrW{1,jj} = 'Thr_W';
% result_EneB{1,jj} = 'Ene_B';
% result_EneW{1,jj} = 'Ene_W';
% 
% 
% outname=sprintf('result-(%g)-(%d-%g-%d)-(%d-%g-%d)-(%d-%d).mat', TR_B, ...
%                         T_B, Rho_B, CW_B(1), ...
%                         T_W, Rho_W, CW_W(1), TUNE_FACTOR, TUNE_FACTOR2);
% NN_W = N_W/2;
% NN_B = N_B/2;
% 
% % load(outname);
% if ~exist(outname, 'file') %
% 
%     Sr_B(1:NN_W+1,1:NN_B+1) = 0;
%     Pe_B(1:NN_W+1,1:NN_B+1) = 0;
%     Sr_W(1:NN_W+1,1:NN_B+1) = 0;
%     En_W(1:NN_W+1,1:NN_B+1) = 0;
% 
% %     for traffic_rate_bmac = [0.2 0.5 1 2] % X = X packets/s
%     for n_B = 0:30 % NN_B % 
%         for n_W = 0:30 % NN_W % 
%             if  n_W~=0 || n_B~=0
%                 out = model(-1, n_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                                         n_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1));
%             end
%         end
% %                                         P_free(n) = out(1);
% %                                         Pe(n) = out(2);
% %                                         Service2(n) = out(3);
% %                                         Service(n) = out(4);
% %         fprintf('beta %i: %g\n', n, nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n));
% %         SS = SS + nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n)*S(n);
%     end
% %                                             Pe_B(NN_W+1,NN_B+1)
% 
% %             return;
% 
%     save(outname, 'Sr_B', 'Pe_B', 'Sr_W', 'En_W', 'P_Wa', 'P_Wd');
% 
% 
% else
%     fprintf('found the file\n');
%     load(outname);
% end  
% 
% % return;
% 
% 
% for N_B = 10:10:30 % 30 % N_B senders
%     
%     jj = jj + 1;
%     result_ThrB{1,jj} = sprintf('N_B=%d mod', N_B);
%     result_ThrW{1,jj} = sprintf('N_B=%d mod', N_B);
%     result_EneB{1,jj} = sprintf('N_B=%d mod', N_B);
%     result_EneW{1,jj} = sprintf('N_B=%d mod', N_B);
% 
%     ii = 1;
%     
%     for N_W = 5:5:30 % [30] % N_W senders
%         NN_W = N_W;
%         NN_B = N_B;
% 
%         fprintf('=============================================\n');
%         fprintf('# LPL Nodes: %d (%d senders)\n', N_B*2, N_B);
%         fprintf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
%         fprintf('Traffic: %g p/s\n', TR_B);
%         fprintf('# PSM Nodes: %d (%d senders)\n', N_W*2, N_W);
%         fprintf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
%         fprintf('Traffic: %g p/s\n', TR_W);
%         fprintf('=============================================\n');
%         
%         
%         [SS_B, SS_W, Q0_B, en_W] = compute(NN_B, TR_B, T_B, NN_W, TR_W, T_W);
%         
%         ii = ii + 1;
% 
% 
%         result_ThrB{ii,jj} = 1-SS_B;
%         result_ThrW{ii,jj} = SS_W;
%         result_EneB{ii,jj} = 1-Q0_B;
%         result_EneW{ii,jj} = en_W;
%         
%         result_ThrB{ii,1} = sprintf('N_W=%d', N_W);
%         result_ThrW{ii,1} = sprintf('N_W=%d', N_W);
%         result_EneB{ii,1} = sprintf('N_W=%d', N_W);
%         result_EneW{ii,1} = sprintf('N_W=%d', N_W);
% 
% 
% 
%     end
% end
% 
% return;



clear all;
clear all;
clear all;



global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
global NN_B NN_W;
global TR_B TR_W;
global T_B T_W;
global q0b q0w;
global CASE;
global TUNE_FACTOR TUNE_FACTOR2;


%defaults:
L_W = 150; %Packet size in DATA us
Ackto_W = 10;

cw_W = 16;
CW_W(1:2) = [cw_W 16];%minimal CW size

N_W = 40; % N_W/2 senders
T_W = 100*1e3;% ms
Rho_W = 0.1;
TR_W = 5;

L_B = 1020;%3000;%us
Ackto_B = 300;%256; %us

cw_B = 70;
CW_B(1:2) = [cw_B 70];

N_B = 40; % N_W/2 senders
T_B = 400*1e3;% ms
Rho_B = 0.05;
TR_B = 1;

TUNE_FACTOR = 800;%0;%
TUNE_FACTOR2 = 250;%0;%

%Figure type 2 (vary # nodes)
% x-axis: # PSM (5 10 15 20)
% each line: # LPL (10 20 30)
% y-axis: Throughput of PSM
% y-axis: Energy consumption of PSM
% y-axis: Throughput of LPL
% y-axis: Energy consumption of LPL
ii = 1;
jj = 1;

result_ThrB2{1,jj} = 'Thr_B';
result_ThrW2{1,jj} = 'Thr_W';
result_EneB2{1,jj} = 'Ene_B';
result_EneW2{1,jj} = 'Ene_W';



NN_W = N_W/2;
NN_B = N_B/2;

NN_W_g = 40/2;
NN_B_g = 40/2;


for TR_B = 1%[0.5 1 3] % [2 2.1 2.2 2.3 2.4 2.5] %
    
    outname=sprintf('result-(%g)-(%d-%g-%d)-(%d-%g-%d)-(%d-%d).mat', TR_B, ...
                            T_B, Rho_B, CW_B(1), ...
                            T_W, Rho_W, CW_W(1), TUNE_FACTOR, TUNE_FACTOR2);

%     load(outname);
    
    if ~exist(outname, 'file') %1%

        Sr_B(1:NN_W_g+1,1:NN_B_g+1) = 0;
        Pe_B(1:NN_W_g+1,1:NN_B_g+1) = 0;
        Sr_W(1:NN_W_g+1,1:NN_B_g+1) = 0;
        En_W(1:NN_W_g+1,1:NN_B_g+1) = 0;

    %     for traffic_rate_bmac = [0.2 0.5 1 2] % X = X packets/s
        for n_B = 0:NN_B_g % NN_B % 
            for n_W = 0:NN_W_g % NN_W % 
                if  n_W~=0 || n_B~=0
                    out = model(-1, n_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                            n_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1));
                end
            end
    %                                         P_free(n) = out(1);
    %                                         Pe(n) = out(2);
    %                                         Service2(n) = out(3);
    %                                         Service(n) = out(4);
    %         fprintf('beta %i: %g\n', n, nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n));
    %         SS = SS + nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n)*S(n);
        end
    %                                             Pe_B(NN_W+1,NN_B+1)

%         return;

        save(outname, 'Sr_B', 'Pe_B', 'Sr_W', 'En_W', 'P_Wa', 'P_Wd');


    else
        fprintf('found the file\n');
        load(outname);
    end  

% return;    
    
%     continue;
    
    
    jj = jj + 1;
    result_ThrB2{1,jj} = sprintf('TR_B=%g mod', TR_B);
    result_ThrW2{1,jj} = sprintf('TR_B=%g mod', TR_B);
    result_EneB2{1,jj} = sprintf('TR_B=%g mod', TR_B);
    result_EneW2{1,jj} = sprintf('TR_B=%g mod', TR_B);

    ii = 1;
    
    for TR_W = 5%1:2:11 % [1 11]% [1 3 5]%

        fprintf('=============================================\n');
        fprintf('# LPL Nodes: %d (%d senders)\n', N_B*2, N_B);
        fprintf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
        fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
        fprintf('Traffic: %g p/s\n', TR_B);
        fprintf('# PSM Nodes: %d (%d senders)\n', N_W*2, N_W);
        fprintf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
        fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
        fprintf('Traffic: %g p/s\n', TR_W);
        fprintf('=============================================\n');
        
        [SS_B, SS_W, Q0_B, en_W] = compute(NN_B, TR_B, T_B, NN_W, TR_W, T_W);

        ii = ii + 1;

        result_ThrB2{ii,jj} = 1-SS_B;
        result_ThrW2{ii,jj} = SS_W;
        result_EneB2{ii,jj} = 1-Q0_B;
        result_EneW2{ii,jj} = en_W;

        
        result_ThrB2{ii,1} = sprintf('TR_W=%g', TR_W);
        result_ThrW2{ii,1} = sprintf('TR_W=%g', TR_W);
        result_EneB2{ii,1} = sprintf('TR_W=%g', TR_W);
        result_EneW2{ii,1} = sprintf('TR_W=%g', TR_W);

        y = (NN_B*(1-Q0_B) + NN_B*Rho_B + N_W*en_W)/(N_W+N_B)


    end
end

return;

clear all;
clear all;
clear all;



global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
global NN_B NN_W;
global TR_B TR_W;
% global T_B T_W;
global q0b q0w;
global CASE;
global TUNE_FACTOR TUNE_FACTOR2;



%defaults:
L_W = 150; %Packet size in DATA us
Ackto_W = 10;

cw_W = 16;
CW_W(1:2) = [cw_W 16];%minimal CW size

N_W = 40; % N_W/2 senders
T_W = 100*1e3;% ms
Rho_W = 0.1;
TR_W = 5;

L_B = 1020;%3000;%us
Ackto_B = 300;%256; %us

cw_B = 70;
CW_B(1:2) = [cw_B 70];

N_B = 40; % N_W/2 senders
T_B = 400*1e3;% ms
Rho_B = 0.05;
TR_B = 1;

TUNE_FACTOR = 800;%0;%
TUNE_FACTOR2 = 250;%0;%

%Figure type 3 (vary # nodes)
% x-axis: # PSM (5 10 15 20)
% each line: # LPL (10 20 30)
% y-axis: Throughput of PSM
% y-axis: Energy consumption of PSM
% y-axis: Throughput of LPL
% y-axis: Energy consumption of LPL
% ii = 1;
% jj = 1;

result_ThrB2{1,1} = 'Thr_B';
result_ThrW2{1,1} = 'Thr_W';
result_EneB2{1,1} = 'Ene_B';
result_EneW2{1,1} = 'Ene_W';



NN_W = N_W/2;
NN_B = N_B/2;

NN_W_g = 60/2;
NN_B_g = 60/2;


fprintf('=============================================\n');
fprintf('# LPL Nodes: %d (%d senders)\n', N_B*2, N_B);
fprintf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
% fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
fprintf('Traffic: %g p/s\n', TR_B);
fprintf('# PSM Nodes: %d (%d senders)\n', N_W*2, N_W);
fprintf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
% fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
fprintf('Traffic: %g p/s\n', TR_W);
fprintf('=============================================\n');


T_B=[200 400 800];%200;%800;%
Rho_B=[0.1 0.04 0.0125];%0.1;%0.0125;%
T_W=[50  50  100  100 200 200];%200;% 200;%200;  
Rho_W=[0.1 0.2 0.05 0.1 0.025 0.05];%0.025;%0.025;%0.1;%   
T_B = T_B.*1e3;
T_W = T_W.*1e3;

% T_W = 100*1e3;% ms
% Rho_W = 0.1;
% T_B = 400*1e3;% ms
% Rho_B = 0.05;






for jj=1:size(T_B,2)
    
    result_ThrB2{1,jj+1} = sprintf('T_B=%g, CR_B=%g mod', T_B(jj)/1e3, Rho_B(jj));
    result_ThrW2{1,jj+1} = sprintf('T_B=%g, CR_B=%g mod', T_B(jj)/1e3, Rho_B(jj));
    result_EneB2{1,jj+1} = sprintf('T_B=%g, CR_B=%g mod', T_B(jj)/1e3, Rho_B(jj));
    result_EneW2{1,jj+1} = sprintf('T_B=%g, CR_B=%g mod', T_B(jj)/1e3, Rho_B(jj));
    
    
    for ii=1:size(T_W,2)
        
        fprintf('LPL cycle: %d ms, duty-cycle ratio: %g\n', T_B(jj), Rho_B(jj));
        fprintf('PSM cycle: %d ms, duty-cycle ratio: %g\n', T_W(ii), Rho_W(ii));

%         outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat\n', sim_time/1000, ...
%                                     N_B, L_B, Ackto_B, T_B(jj), Rho_B(jj), CW_B(1), TR_B, ...
%                                     N_W, L_W, Ackto_W, T_W(ii), Rho_W(ii), CW_W(1), TR_W);
                                
        outname=sprintf('result-(%g)-(%d-%g-%d)-(%d-%g-%d)-(%d-%d).mat', TR_B, ...
                                T_B(jj), Rho_B(jj), CW_B(1), ...
                                T_W(ii), Rho_W(ii), CW_W(1), TUNE_FACTOR, TUNE_FACTOR2);
                                
%         disp(outname);
%         continue;
                              
        if 1%~exist(outname, 'file') %

            Sr_B(1:NN_W_g+1,1:NN_B_g+1) = 0;
            Pe_B(1:NN_W_g+1,1:NN_B_g+1) = 0;
            Sr_W(1:NN_W_g+1,1:NN_B_g+1) = 0;
            En_W(1:NN_W_g+1,1:NN_B_g+1) = 0;

        %     for traffic_rate_bmac = [0.2 0.5 1 2] % X = X packets/s
            for n_B = 0:NN_B_g % NN_B % 
                for n_W = 0:NN_W_g % NN_W % 
                    if  n_W~=0 || n_B~=0
                        out = model(-1, n_B, L_B, Ackto_B, T_B(jj), Rho_B(jj), CW_B(1), TR_B, ...
                                                n_W, L_W, Ackto_W, T_W(ii), Rho_W(ii), CW_W(1));
                    end
                end
%                                         P_free(n) = out(1);
%                                         Pe(n) = out(2);
%                                         Service2(n) = out(3);
%                                         Service(n) = out(4);
        %         fprintf('beta %i: %g\n', n, nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n));
        %         SS = SS + nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n)*S(n);
            end
%                                             Pe_B(NN_W+1,NN_B+1)

            return;

            save(outname, 'Sr_B', 'Pe_B', 'Sr_W', 'En_W', 'P_Wa', 'P_Wd');


        else
%             fprintf('found the file\n');
            load(outname);
        end        
%         disp(outname);
%         continue;
%         load(outname, 'node_struct')
%         [Pe_B, Th_W, ene_B, ene_W] = compute(sim_time, N_W, N_B, node_struct);
%         NN_B
%         NN_W

        [SS_B, SS_W, Q0_B, en_W] = compute(NN_B, TR_B, T_B(jj), NN_W, TR_W, T_W(ii));


%         SS_B = 0.1;
%         SS_W = 0.2;
%         Q0_B = 0.3;
%         en_W = 0.4;

        result_ThrB2{ii+1,jj+1} = 1-SS_B;
        result_ThrW2{ii+1,jj+1} = SS_W;
        result_EneB2{ii+1,jj+1} = 1-Q0_B;
        result_EneW2{ii+1,jj+1} = en_W;
        
        fprintf('result_ThrB2{%d,%d}: %g\n', ii+1, jj+1, result_ThrB2{ii+1,jj+1});
        fprintf('result_ThrW2{%d,%d}: %g\n', ii+1, jj+1, result_ThrW2{ii+1,jj+1});
        fprintf('result_EneB2{%d,%d}: %g\n', ii+1, jj+1, result_EneB2{ii+1,jj+1});
        fprintf('result_EneW2{%d,%d}: %g\n', ii+1, jj+1, result_EneW2{ii+1,jj+1});

        result_ThrB2{ii+1,1} = sprintf('T_W=%g, CR_W=%g mod', T_W(ii)/1e3, Rho_W(ii));
        result_ThrW2{ii+1,1} = sprintf('T_W=%g, CR_W=%g mod', T_W(ii)/1e3, Rho_W(ii));
        result_EneB2{ii+1,1} = sprintf('T_W=%g, CR_W=%g mod', T_W(ii)/1e3, Rho_W(ii));
        result_EneW2{ii+1,1} = sprintf('T_W=%g, CR_W=%g mod', T_W(ii)/1e3, Rho_W(ii));
%         return;
    end
    
end




% return;
% 
% clear all;
% clear all;
% clear all;
% 
% global aMinBE LDATA;
% global showdetailed_nt;
% global showdetailed_st;
% global showdetailed_ct;
% showdetailed_nt = 0;
% showdetailed_st = 0;
% showdetailed_ct = 0;
% global savefigure;
% savefigure = 0;
% 
% if showdetailed_nt
%     global afig singlefig;
%     singlefig = 1;
%     if singlefig
%         afig = figure;
%     end
% end
% if showdetailed_st
%     global sfig;
%     sfig = figure;
% end
% if showdetailed_ct
%     global cfig;
%     cfig = figure;
% end
% 
% 
% cwfig = figure;
% 
% % global Stage_cur Counter_cur Actnode_cur;
% global CW0_stat CW1_stat CW2_stat CW3_stat CW4_stat CW5_stat CW_stat;
% col1=hsv(100);
% global Actnode_dist;
% 
% kk=0;
% 
% for aMinBE=16:4:16%[20 24 28 36 40 44 48 52 56 60]%
% 
% 
% zw_time(1:6)=0;
% round = 1000; %set round to simulate
% % T0 = 1/1;%s
% % T = T0*100000;
% BI = 500;
% T = round*BI;
% T0 = T/100000;
% global col;
% col=hsv(round);
% 
% 
% PAYLOAD = 500;%1200;%200;%
% 
% for atim=200:200:200 %for Daren Cline's test, get rid of ATIM
%     data=BI-atim;
%     for n= 10:5:10%[38 40 43 45 48] %[20 22 25 28 30 32 35]%[50 53 55] %[10 13 15 17]% 
%         zw_time=clock();
%         outname=sprintf('mod(%d,%d,%d)-%d-%d-%d-%d-%d-%g.txt', ...
%             n, atim, data, zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
%         diary(outname);
%         diary on;
%         for tr=50:2.5:50
%             diary on;
%             [stat1,stat2,stat3, pp, round] = lpl(T, n, PAYLOAD, tr, atim, data);
%             diary on;
%             fprintf('\n\n');
%             diary off;
%             kk = kk + 1;
%         end
% 
%         figure;
%         plot(stat3/(round), 'DisplayName', 'A(t)');
%         xlabel('Time');ylabel('Average A(t)');
%         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
% %         hold on;
%         if savefigure
%             outname=sprintf('N(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
%             savefig(outname);
%         end
% %         % 
% %         figure;
% %         plot(stat1/(round), 'DisplayName', 'S(t)');
% %         xlabel('Time');ylabel('Average S(t)');
% %         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
% % %         hold on;
% %         if savefigure
% %             outname=sprintf('S(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
% %             savefig(outname);
% %         end
% %         % 
% %         figure;
% %         plot(stat2/(round), 'DisplayName', 'C(t)');
% %         xlabel('Time');ylabel('Average C(t)');
% %         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
% % %         hold on;
% %         if savefigure
% %             outname=sprintf('C(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
% %             savefig(outname);
% %         end
% 
%         outname=sprintf('stat(%d-%d-%d-%d-%d).mat', round, n, aMinBE, BI, LDATA);
% 
%         save(outname, 'stat1', 'stat2', 'stat3', 'round', 'CW_stat', 'Actnode_dist');
% 
%         
%         figure;
%         plot(log2(CW_stat/(round)), 'DisplayName', 'CW(t)');
%         hold on;
%         
% %         figure(cwfig);
% %         plot(log2(CW_stat/(round)), 'DisplayName', sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round), 'Color',col1(kk,:));
% %         xlabel('Time');ylabel('Average CW(t)');
% % %         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, 2^aMinBE, BI, LDATA, round));
% %         hold on;
% %         
% %         if savefigure
% %             outname=sprintf('CW(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
% %             savefig(outname);
% %         end
% 
% %         close('all'); %close all figures
%     end
% end
% 
% end
% 
% figure;
% %  colormap('hot');   % set colormap
% imagesc(Actnode_dist/round);        % draw image and scale colormap to values range
% colorbar;          % show color scale
% 
% 
% % if round==1
% % 
% %    
% %     for t=1:BI
% %         CW0_stat(t) = size(find(Stage_cur(:,t)==0),1)/Actnode_cur(t);
% %         CW1_stat(t) = size(find(Stage_cur(:,t)==1),1)/Actnode_cur(t);
% %         CW2_stat(t) = size(find(Stage_cur(:,t)==2),1)/Actnode_cur(t);
% %         CW3_stat(t) = size(find(Stage_cur(:,t)==3),1)/Actnode_cur(t);
% %         CW4_stat(t) = size(find(Stage_cur(:,t)==4),1)/Actnode_cur(t);
% %     end
% %     
% % end


