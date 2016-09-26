% y = randn(1,1000);
% h = plot(y(1:200));
% ylim([-5 5]);
% n = 2;
% while n+199 <= numel(y)
%    newy = y(n+(1:199));
%    disp(newy);
%    set(h,'ydata',newy);
%    drawnow;
%    n = n+1;
%    return;
% end
% 
% return;

% T = 100;%s
% T = T*100000;
% n = 30;
% tr = 10;
% PAYLOAD = 1500;
% atim = 100;
% data = 2000 - atim - 1;
% 
% % zw_time=clock();
% % outname=['temp-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
% %     num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
% %     '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];
% % diary(outname);
% % diary on;
% 
% wifi_psm(T, n, PAYLOAD, tr, atim, data);
% % md1(T, n, PAYLOAD, tr, atim, data);
% 
% diary off;
% return;

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



for N_W = 40 %  10:10:60 % N_W/2 senders
    
    L_W = 15+2;% WiFi time slots (15 == 1000 bytes)
    Ackto_W = 4;% WiFi time slots
    for T_W = 200 %50 %  ms
        for Rho_W = [0.025 0.05]%0.2 %
            CW_W(1:2) = [16 16];%minimal CW size
            for TR_W = 5 %1:2:11 % X = X packets/s


                round = 100; %set round to simulate
                global col;
                col=hsv(round);

                % T0 = T/100000;

                sim_time = round*T_W;%10000;%ms


            % for N = 40:20:40
            % %odd id for sender
            % %even id for receiver.
            % %receiver_id = sender_id + 1

                % t = 0;
                L_B = 34; %10; %100;%BMAC time slots (100 == 94 bytes)
                Ackto_B = 10;%256; % BMAC time slots (= 300us)
                % timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot


                % mac_queue(1:N) = 0; %# packets
                % QUEUE_SIZE = 100;

                for T_B = 200 %400:200:400 %ms
                    for Rho_B = 0.1 %[0.05 0.1]%
                        % T_l = T*rho;
                        % T_s = T-T_l;

                        % timeslot = 10;%us
                        % tic;
                        % true_time = sim_time*1e3/timeslot;

                        for cw_B=70 %[70 310]%

                            CW_B(1:2) = [cw_B 70];
                %             CW(1:2) = [310 70];
        %                     jj = jj + 1;
        %                     ii = 0;
                            for TR_B = 1 %[0.5 1 3] % [0.1 0.2] %% X = X packets/s
                                for N_B = 40%[10 20]%10:10:60 %N_B/2 senders
                %odd id for sender
                %even id for receiver.
                %receiver_id = sender_id + 1

                                    zw_printf('=============================================\n');
                                    zw_printf('Simulation time: %g s\n', sim_time/1000);
                                    zw_printf('# LPL Nodes: %d (%d senders)\n', N_B, N_B/2);
                                    zw_printf('CW: %d time slots, Packet size: %d\n', CW_B(1), L_B);
                                    zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
                                    zw_printf('Traffic: %g p/s\n', TR_B);
                                    zw_printf('# PSM Nodes: %d (%d senders)\n', N_W, N_W/2);
                                    zw_printf('CWmin: %d time slots, Packet size: %d\n', CW_W(1), L_W);
                                    zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
                                    zw_printf('Traffic: %g p/s\n', TR_W);
                                    zw_printf('=============================================\n');

                                    outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...% 500, ...%
                                                                N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                                                N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
                                    if 1%~exist(outname, 'file') %
                                        
                                        node_struct = cox(sim_time, N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                                                    N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);

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
                                                                    N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                                                    N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
                                        savefig(outname);
%                                         close(asfig);

                                    else
                                        fprintf('found the file\n');
                                        load(outname, 'node_struct', 'packet_mactime', 'packet_sucinatim');
                                    end
            %                         clear node_struct;
            %                         load(outname, 'node_struct');

                                    firstmactime(N_W/2,1) = sum(packet_mactime(:,1))/round;
                                    zw_printf('First packet Sertime(%d,%d): %g\n\n', N_B/2, N_W/2, sum(packet_mactime(:,1))/round);

                                    Pe_B(N_W/2,N_B/2) = 1-(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B/2)/(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B + sum([node_struct(N_W/2+1:N_W/2+N_B).drop]));
                                    ii = ii + 1;

                                    zw_printf('Pe_B(%d,%d): %g\n\n', N_B/2, N_W/2, Pe_B(N_W/2,N_B/2));
                                    zw_printf('Th_B(%d,%d): %g\n\n', N_B/2, N_W/2, 1-Pe_B(N_W/2,N_B/2));
                                    
%                                     if ismember(N_W, [10:10:40]) && ismember(TR_W, [1:2:7]) && ismember(N_B, [10:10:40]) && ismember(TR_B, [0.1 0.5 1 2])
%                                         sim_time = 5000*T_W;
%                                     else
%                                         sim_time = 1000*T_W;
%                                     end

                                    Th_W(N_W/2,N_B/2) = (sum([node_struct(1:N_W/2).success]))/(N_W/2)/(sim_time/1e3);
                                    zw_printf('Th_W(%d,%d): %g p/s\n\n', N_B/2, N_W/2, Th_W(N_W/2,N_B/2));

            %                         thr(ii,jj) = (sum([node_struct.success]) + N_B/2)/(sum([node_struct.success]) + N_B + sum([node_struct.drop]));
            % 
            %                         zw_printf('Throughput: %g\n\n', thr(ii,jj));

                %                     total_txtime = sum([node_struct.txtime])/(sum([node_struct.count]))
                %                     find([node_struct.txtime]>=40000)
                                    ene_B(N_W/2,N_B/2) = sum([node_struct(N_W/2+1:N_W/2+N_B).awaketime])/(N_B/2)/100/sim_time;
                                    zw_printf('Ene_B(%d,%d): %g\n\n', N_B/2, N_W/2, ene_B(N_W/2,N_B/2));

                                    ene_W(N_W/2,N_B/2) = sum([node_struct(1:N_W/2).awaketime])/(N_W/2)/100/sim_time;
                                    zw_printf('Ene_W(%d,%d): %g\n\n', N_B/2, N_W/2, ene_W(N_W/2,N_B/2));


                                    result(ii,1) = 1-Pe_B(N_W/2,N_B/2);
                                    result(ii,2) = Th_W(N_W/2,N_B/2);
                                    result(ii,3) = ene_B(N_W/2,N_B/2);
                                    result(ii,4) = ene_W(N_W/2,N_B/2);
                                    zzww{ii,1} = sprintf('N_W:%d,TR_W:%d,N_B:%d,TR_B:%g', N_W/2, TR_W, N_B/2, TR_B);

%                                     return;

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
return;






















clear all;
clear all;
clear all;

global aMinBE LDATA;
global showdetailed_nt;
global showdetailed_st;
global showdetailed_ct;
showdetailed_nt = 0;
showdetailed_st = 0;
showdetailed_ct = 0;
global savefigure;
savefigure = 0;

if showdetailed_nt
    global afig singlefig;
    singlefig = 1;
    if singlefig
        afig = figure;
    end
end
if showdetailed_st
    global sfig;
    sfig = figure;
end
if showdetailed_ct
    global cfig;
    cfig = figure;
end


cwfig = figure;

% global Stage_cur Counter_cur Actnode_cur;
global CW0_stat CW1_stat CW2_stat CW3_stat CW4_stat CW5_stat CW_stat;
col1=hsv(100);
global Actnode_dist;

kk=0;

for aMinBE=16:4:16%[20 24 28 36 40 44 48 52 56 60]%


zw_time(1:6)=0;
round = 1000; %set round to simulate
% T0 = 1/1;%s
% T = T0*100000;
BI = 500;
T = round*BI;
T0 = T/100000;
global col;
col=hsv(round);


PAYLOAD = 500;%1200;%200;%

for atim=200:200:200 %for Daren Cline's test, get rid of ATIM
    data=BI-atim;
    for n= 10:5:10%[38 40 43 45 48] %[20 22 25 28 30 32 35]%[50 53 55] %[10 13 15 17]% 
        zw_time=clock();
        outname=sprintf('mod(%d,%d,%d)-%d-%d-%d-%d-%d-%g.txt', ...
            n, atim, data, zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
        diary(outname);
        diary on;
        for tr=50:2.5:50
            diary on;
            [stat1,stat2,stat3, pp, round] = lpl(T, n, PAYLOAD, tr, atim, data);
            diary on;
            fprintf('\n\n');
            diary off;
            kk = kk + 1;
        end

        figure;
        plot(stat3/(round), 'DisplayName', 'A(t)');
        xlabel('Time');ylabel('Average A(t)');
        title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
%         hold on;
        if savefigure
            outname=sprintf('N(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
            savefig(outname);
        end
%         % 
%         figure;
%         plot(stat1/(round), 'DisplayName', 'S(t)');
%         xlabel('Time');ylabel('Average S(t)');
%         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
% %         hold on;
%         if savefigure
%             outname=sprintf('S(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
%             savefig(outname);
%         end
%         % 
%         figure;
%         plot(stat2/(round), 'DisplayName', 'C(t)');
%         xlabel('Time');ylabel('Average C(t)');
%         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round));
% %         hold on;
%         if savefigure
%             outname=sprintf('C(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
%             savefig(outname);
%         end

        outname=sprintf('stat(%d-%d-%d-%d-%d).mat', round, n, aMinBE, BI, LDATA);

        save(outname, 'stat1', 'stat2', 'stat3', 'round', 'CW_stat', 'Actnode_dist');

        
        figure;
        plot(log2(CW_stat/(round)), 'DisplayName', 'CW(t)');
        hold on;
        
%         figure(cwfig);
%         plot(log2(CW_stat/(round)), 'DisplayName', sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, aMinBE, BI, LDATA, round), 'Color',col1(kk,:));
%         xlabel('Time');ylabel('Average CW(t)');
% %         title(sprintf('N=%d, W_0=%d, BI=%d, d=%d, over %d rounds', n, 2^aMinBE, BI, LDATA, round));
%         hold on;
%         
%         if savefigure
%             outname=sprintf('CW(t)-%d-%d-%d-%d-%d-new.fig', round, n, aMinBE, BI, LDATA);
%             savefig(outname);
%         end

%         close('all'); %close all figures
    end
end

end

figure;
%  colormap('hot');   % set colormap
imagesc(Actnode_dist/round);        % draw image and scale colormap to values range
colorbar;          % show color scale


% if round==1
% 
%    
%     for t=1:BI
%         CW0_stat(t) = size(find(Stage_cur(:,t)==0),1)/Actnode_cur(t);
%         CW1_stat(t) = size(find(Stage_cur(:,t)==1),1)/Actnode_cur(t);
%         CW2_stat(t) = size(find(Stage_cur(:,t)==2),1)/Actnode_cur(t);
%         CW3_stat(t) = size(find(Stage_cur(:,t)==3),1)/Actnode_cur(t);
%         CW4_stat(t) = size(find(Stage_cur(:,t)==4),1)/Actnode_cur(t);
%     end
%     
% end


