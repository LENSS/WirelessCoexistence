
clear all;
clear all;
clear all;

global DEBUG_PRT DEBUG_LOG;
DEBUG_PRT = 1;
DEBUG_LOG = 0; % log to file
ii = 0;
jj = 0;

global showdetailed_nt;
showdetailed_nt = 0;

global afig singlefig;
global asfig;

if showdetailed_nt
    singlefig = 1;
    if singlefig
        afig = figure;
    end
end

global packet_mactime;
global packet_sucinatim;

% defaults;
% 
% 
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
% for N_B = 60 %20:20:60 %N_B/2 senders
%     jj = jj + 1;
%     result_ThrB{1,jj} = sprintf('N_B=%d sim', N_B/2);
%     result_ThrW{1,jj} = sprintf('N_B=%d sim', N_B/2);
%     result_EneB{1,jj} = sprintf('N_B=%d sim', N_B/2);
%     result_EneW{1,jj} = sprintf('N_B=%d sim', N_B/2);
% 
%     ii = 1;
%     for N_W = 60 %10:10:60 %N_W/2 senders
%  
%         fprintf('=============================================\n');
%         fprintf('Simulation time: %g s\n', sim_time/1000);
%         fprintf('# LPL Nodes: %d (%d senders)\n', N_B, N_B/2);
%         fprintf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
%         fprintf('Traffic: %g p/s\n', TR_B);
%         fprintf('# PSM Nodes: %d (%d senders)\n', N_W, N_W/2);
%         fprintf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
%         fprintf('Traffic: %g p/s\n', TR_W);
%         fprintf('=============================================\n');
% 
%         sim_time = 5000*T_W;
% 
%         outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...
%                                     N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                                     N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
% 
%         load(outname, 'node_struct')
%         
%         if ismember(N_W, [10:10:40]) && ismember(TR_W, [1:2:7]) && ismember(N_B, [10:10:40]) && ismember(TR_B, [0.1 0.5 1 2])
%             sim_time = 5000*T_W;
%         else
%             sim_time = 1000*T_W;
%         end
% 
%         Pe_B(N_W/2,N_B/2) = 1-(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B/2)/(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B + sum([node_struct(N_W/2+1:N_W/2+N_B).drop]));
%         ii = ii + 1;
% 
%         fprintf('Pe_B(%d,%d): %g\n\n', N_B/2, N_W/2, Pe_B(N_W/2,N_B/2));
%         fprintf('Th_B(%d,%d): %g\n\n', N_B/2, N_W/2, 1-Pe_B(N_W/2,N_B/2));
% 
%         Th_W(N_W/2,N_B/2) = (sum([node_struct(1:N_W/2).success]))/(N_W/2)/(sim_time/1e3);
%         fprintf('Th_W(%d,%d): %g p/s\n\n', N_B/2, N_W/2, Th_W(N_W/2,N_B/2));
% 
%         ene_B(N_W/2,N_B/2) = sum([node_struct(N_W/2+1:N_W/2+N_B).awaketime])/(N_B/2)/100/sim_time;
%         fprintf('Ene_B(%d,%d): %g\n\n', N_B/2, N_W/2, ene_B(N_W/2,N_B/2));
% 
%         ene_W(N_W/2,N_B/2) = sum([node_struct(1:N_W/2).awaketime])/(N_W/2)/100/sim_time;
%         fprintf('Ene_W(%d,%d): %g\n\n', N_B/2, N_W/2, ene_W(N_W/2,N_B/2));
% 
%         result_ThrB{ii,jj} = 1-Pe_B(N_W/2,N_B/2);
%         result_ThrW{ii,jj} = Th_W(N_W/2,N_B/2);
%         result_EneB{ii,jj} = ene_B(N_W/2,N_B/2);
%         result_EneW{ii,jj} = ene_W(N_W/2,N_B/2);
%         
%         result_ThrB{ii,1} = sprintf('N_W=%d', N_W/2);
%         result_ThrW{ii,1} = sprintf('N_W=%d', N_W/2);
%         result_EneB{ii,1} = sprintf('N_W=%d', N_W/2);
%         result_EneW{ii,1} = sprintf('N_W=%d', N_W/2);
% 
% 
% 
%     end
% end
% 
% 
% return;

% clear all;


% defaults;
% 
%     
%     
% %Figure type 2 (vary data arrival rate)
% % x-axis: TR PSM (1 3 5 7 9 11)
% % each line: # LPL (0.5 1 3)
% % y-axis: Throughput of PSM
% % y-axis: Energy consumption of PSM
% % y-axis: Throughput of LPL
% % y-axis: Energy consumption of LPL
% 
% ii = 1;
% jj = 1;
% 
% result_ThrB2{1,jj} = 'Thr_B';
% result_ThrW2{1,jj} = 'Thr_W';
% result_EneB2{1,jj} = 'Ene_B';
% result_EneW2{1,jj} = 'Ene_W';
% 
% for TR_B = [0.5 1 3] %N_W/2 senders
%     jj = jj + 1;
%     result_ThrB2{1,jj} = sprintf('TR_B=%g sim', TR_B);
%     result_ThrW2{1,jj} = sprintf('TR_B=%g sim', TR_B);
%     result_EneB2{1,jj} = sprintf('TR_B=%g sim', TR_B);
%     result_EneW2{1,jj} = sprintf('TR_B=%g sim', TR_B);
% 
%     ii = 1;
%     for TR_W = 1:2:11 %N_B/2 senders
%  
%         fprintf('=============================================\n');
%         fprintf('Simulation time: %g s\n', sim_time/1000);
%         fprintf('# LPL Nodes: %d (%d senders)\n', N_B, N_B/2);
%         fprintf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
%         fprintf('Traffic: %g p/s\n', TR_B);
%         fprintf('# PSM Nodes: %d (%d senders)\n', N_W, N_W/2);
%         fprintf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
%         fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
%         fprintf('Traffic: %g p/s\n', TR_W);
%         fprintf('=============================================\n');
% 
%         sim_time = 10000*T_W;
% 
%         outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...
%                                     N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                                     N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
% 
%         load(outname, 'node_struct')
%         
% %         if ismember(N_W, [10:10:40]) && ismember(TR_W, [1:2:7]) && ismember(N_B, [10:10:40]) && ismember(TR_B, [0.1 0.5 1 2])
% %             sim_time = 5000*T_W;
% %         else
% %             sim_time = 1000*T_W;
% %         end
%         
%         Pe_B(N_W/2,N_B/2) = 1-(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B/2)/(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B + sum([node_struct(N_W/2+1:N_W/2+N_B).drop]));
%         ii = ii + 1;
% 
%         fprintf('Pe_B(%d,%d): %g\n\n', N_B/2, N_W/2, Pe_B(N_W/2,N_B/2));
%         fprintf('Th_B(%d,%d): %g\n\n', N_B/2, N_W/2, 1-Pe_B(N_W/2,N_B/2));
% 
%         Th_W(N_W/2,N_B/2) = (sum([node_struct(1:N_W/2).success]))/(N_W/2)/(sim_time/1e3);
%         fprintf('Th_W(%d,%d): %g p/s\n\n', N_B/2, N_W/2, Th_W(N_W/2,N_B/2));
% 
%         ene_B(N_W/2,N_B/2) = sum([node_struct(N_W/2+1:N_W/2+N_B).awaketime])/(N_B/2)/100/sim_time;
%         fprintf('Ene_B(%d,%d): %g\n\n', N_B/2, N_W/2, ene_B(N_W/2,N_B/2));
% 
%         ene_W(N_W/2,N_B/2) = sum([node_struct(1:N_W/2).awaketime])/(N_W/2)/100/sim_time;
%         fprintf('Ene_W(%d,%d): %g\n\n', N_B/2, N_W/2, ene_W(N_W/2,N_B/2));
% 
%         result_ThrB2{ii,jj} = 1-Pe_B(N_W/2,N_B/2);
%         result_ThrW2{ii,jj} = Th_W(N_W/2,N_B/2);
%         result_EneB2{ii,jj} = ene_B(N_W/2,N_B/2);
%         result_EneW2{ii,jj} = ene_W(N_W/2,N_B/2);
%         
%         result_ThrB2{ii,1} = sprintf('TR_W=%g', TR_W);
%         result_ThrW2{ii,1} = sprintf('TR_W=%g', TR_W);
%         result_EneB2{ii,1} = sprintf('TR_W=%g', TR_W);
%         result_EneW2{ii,1} = sprintf('TR_W=%g', TR_W);
% 
% 
% 
%     end
% end



% defaults;

    
    
%Figure type 3 (vary duty cycle ratio and length)
% x-axis: CR PSM ()
% each line: # LPL ()
% y-axis: Throughput of PSM
% y-axis: Energy consumption of PSM
% y-axis: Throughput of LPL
% y-axis: Energy consumption of LPL


result_ThrB2{1,1} = 'Thr_B';
result_ThrW2{1,1} = 'Thr_W';
result_EneB2{1,1} = 'Ene_B';
result_EneW2{1,1} = 'Ene_W';



sim_time = 10*1e3; %500s
N_W = 40;
L_W = 15+2;% WiFi time slots (15 == 1000 bytes)
Ackto_W = 4;% WiFi time slots
CW_W(1:2) = [16 16];%minimal CW size
TR_W = 8;
N_B = 40;
L_B = 34; %10; %100;%BMAC time slots (100 == 94 bytes)
Ackto_B = 10;%256; % BMAC time slots (= 300us)
CW_B(1:2) = [70 70];
TR_B = 4;

zw_printf('=============================================\n');
zw_printf('Simulation time: %g s\n', sim_time/1000);
zw_printf('# LPL Nodes: %d (%d senders)\n', N_B, N_B/2);
zw_printf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
% zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
zw_printf('Traffic: %g p/s\n', TR_B);
zw_printf('# PSM Nodes: %d (%d senders)\n', N_W, N_W/2);
zw_printf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
% zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
zw_printf('Traffic: %g p/s\n', TR_W);
zw_printf('=============================================\n');

T_B=400;%[200 400 800];%
Rho_B=0.04;%[0.1 0.04 0.0125];%
T_W=200;%[50  50  100  100 200 200];%200;  200;% 
Rho_W=0.025;%[0.1 0.2 0.05 0.1 0.025 0.05];%0.1;%   0.025;%


% T_B = [200 400 800];
% Rho_B=[0.1 0.04 0.0125];
% T_W = [50 50 100 100 200 200];%200;   
% Rho_W=[0.1 0.2 0.05 0.1 0.025 0.05];%0.05;

for jj=1:size(T_B,2)
    
    result_ThrB2{1,jj+1} = sprintf('T_B=%g, CR_B=%g sim', T_B(jj), Rho_B(jj));
    result_ThrW2{1,jj+1} = sprintf('T_B=%g, CR_B=%g sim', T_B(jj), Rho_B(jj));
    result_EneB2{1,jj+1} = sprintf('T_B=%g, CR_B=%g sim', T_B(jj), Rho_B(jj));
    result_EneW2{1,jj+1} = sprintf('T_B=%g, CR_B=%g sim', T_B(jj), Rho_B(jj));
    
    
    for ii=1:size(T_W,2)
        
        zw_printf('LPL cycle: %d ms, duty-cycle ratio: %g\n', T_B(jj), Rho_B(jj));
        zw_printf('PSM cycle: %d ms, duty-cycle ratio: %g\n', T_W(ii), Rho_W(ii));

        outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...
                                    N_B, L_B, Ackto_B, T_B(jj), Rho_B(jj), CW_B(1), TR_B, ...
                                    N_W, L_W, Ackto_W, T_W(ii), Rho_W(ii), CW_W(1), TR_W);
                                
                                
                              
        if 1%~exist(outname, 'file') %

            node_struct = cox(sim_time, N_B, L_B, Ackto_B, T_B(jj), Rho_B(jj), CW_B(1), TR_B, ...
                                    N_W, L_W, Ackto_W, T_W(ii), Rho_W(ii), CW_W(1), TR_W);

            save(outname, 'node_struct', 'packet_mactime', 'packet_sucinatim');

            if showdetailed_nt

                figure(afig);
                outname=sprintf('fig\\afig(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).fig', sim_time/1000, ...
                                            N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                            N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
                savefig(outname);
            end
            figure(asfig);
            outname=sprintf('fig\\asfig(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).fig', sim_time/1000, ...
                                        N_B, L_B, Ackto_B, T_B(jj), Rho_B(jj), CW_B(1), TR_B, ...
                                    N_W, L_W, Ackto_W, T_W(ii), Rho_W(ii), CW_W(1), TR_W);
            savefig(outname);
%                                         close(asfig);

        else
            fprintf('found the file\n');
            load(outname, 'node_struct');
        end        
%         disp(outname);
%         continue;
%         load(outname, 'node_struct')
%         [Pe_B, Th_W, ene_B, ene_W] = compute(sim_time, N_W, N_B, node_struct);

        Pe_B = 1-(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B/2)/(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B + sum([node_struct(N_W/2+1:N_W/2+N_B).drop]));

        fprintf('Pe_B(%d,%d): %g\n\n', N_B/2, N_W/2, Pe_B);
        fprintf('Th_B(%d,%d): %g\n\n', N_B/2, N_W/2, 1-Pe_B);

        Th_W = (sum([node_struct(1:N_W/2).success]))/(N_W/2)/(sim_time/1e3);
        fprintf('Th_W(%d,%d): %g p/s\n\n', N_B/2, N_W/2, Th_W);

        ene_B = sum([node_struct(N_W/2+1:N_W/2+N_B).awaketime])/(N_B/2)/100/sim_time;
        fprintf('Ene_B(%d,%d): %g\n\n', N_B/2, N_W/2, ene_B);

        ene_W = sum([node_struct(1:N_W/2).awaketime])/(N_W/2)/100/sim_time;
        fprintf('Ene_W(%d,%d): %g\n\n', N_B/2, N_W/2, ene_W);

%         Pe_B = 0.1;
%         Th_W = 0.2;
%         ene_B = 0.3;
%         ene_W = 0.4;

        result_ThrB2{ii+1,jj+1} = 1-Pe_B;
        result_ThrW2{ii+1,jj+1} = Th_W;
        result_EneB2{ii+1,jj+1} = ene_B;
        result_EneW2{ii+1,jj+1} = ene_W;

        result_ThrB2{ii+1,1} = sprintf('T_W=%g, CR_W=%g sim', T_W(ii), Rho_W(ii));
        result_ThrW2{ii+1,1} = sprintf('T_W=%g, CR_W=%g sim', T_W(ii), Rho_W(ii));
        result_EneB2{ii+1,1} = sprintf('T_W=%g, CR_W=%g sim', T_W(ii), Rho_W(ii));
        result_EneW2{ii+1,1} = sprintf('T_W=%g, CR_W=%g sim', T_W(ii), Rho_W(ii));
        

    end
    
end






