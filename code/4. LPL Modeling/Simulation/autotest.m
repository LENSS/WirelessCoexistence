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

% for N = 40:20:40
% %odd id for sender
% %even id for receiver.
% %receiver_id = sender_id + 1

    % t = 0;
    % L_bmac = 100;%bytes
    Ack_timeout_bamc = 10;%256; %time slots (= 30us)
    % timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot


    % mac_queue(1:N) = 0; %# packets
    % QUEUE_SIZE = 100;

    for T = 400:200:400 %ms
        for rho = 0.1 %[0.05 0.1]%
            % T_l = T*rho;
            % T_s = T-T_l;

            % timeslot = 10;%us
            % tic;
            sim_time = 100000;%10000;%ms
            % true_time = sim_time*1e3/timeslot;
            
            for cw=70 %[70 310]%
                
                CW(1:2) = [cw 70];
    %             CW(1:2) = [310 70];
                jj = jj + 1;
                ii = 0;
                for traffic_rate_bmac = 20%0.5 %[0.1 0.2 0.5 1 2] %[0.1 0.2] %% X = X packets/s
                    for N = 40:20:40
    %odd id for sender
    %even id for receiver.
    %receiver_id = sender_id + 1

                        zw_printf('=============================================\n');
                        zw_printf('Simulation time: %g s\n', sim_time/1000);
                        zw_printf('# Nodes: %d (%d senders)\n', N, N/2);
                        zw_printf('CW: %d time slots\n', CW(1));
                        zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T, rho);
                        zw_printf('Traffic: %g p/s\n', traffic_rate_bmac);
                        zw_printf('=============================================\n');

                        node_struct = lpl(sim_time, N, Ack_timeout_bamc, T, rho, CW(1), traffic_rate_bmac);

                        outname=sprintf('stat(%d-%d-%d-%g-%d-%g).mat', sim_time/1000, N, T, rho, CW(1), traffic_rate_bmac);
                        clear node_struct;
                        load(outname, 'node_struct');

                        ii = ii + 1;
                        thr(ii,jj) = (sum([node_struct.success]) + N/2)/(sum([node_struct.success]) + N + sum([node_struct.drop]));

                        zw_printf('Pe: %g\n\n', 1-thr(ii,jj));
                        zw_printf('Throughput: %g\n\n', thr(ii,jj));

    %                     total_txtime = sum([node_struct.txtime])/(sum([node_struct.count]))
    %                     find([node_struct.txtime]>=40000)
                        ene(ii,jj) = sum([node_struct.awaketime])/(N/2)/100/sim_time;
                        zw_printf('Energy: %g\n\n', ene(ii,jj));

        %                 return;

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


