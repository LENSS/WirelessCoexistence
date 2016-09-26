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

globals;
global Service;

ii = 0;
jj = 0;



% for N = 20:10:20 %# senders
% 
% Pe(1:N) = 0;
% P_free(1:N) = 0;
% Service(1:N) = 0;
% S(1:N) = 0;


    % t = 0;
    % L_bmac = 100;%bytes
    Ack_timeout_bamc = 300;%256; %us
    % timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot


    % mac_queue(1:N) = 0; %# packets
    % QUEUE_SIZE = 100;

    for T = 400000:200000:400000%us
        for rho = 0.05 %[0.05 0.1] %
            % T_l = T*rho;
            % T_s = T-T_l;

            % timeslot = 10;%us
            % tic;
            sim_time = 2000000;%ms
            % true_time = sim_time*1e3/timeslot;
            
            for cw = 70 %[70 310] %

            %     CW(1:2) = [70 70];
                CW(1:2) = [cw 70];
                
                jj = jj + 1;
                ii = 0;

                for traffic_rate_bmac = 2 %[0.1 0.2 0.5 1 2] %[0.1 0.2] %
                    for N = 40%20:10:20 % <===================== # senders

                    Pe(1:N) = 0;
                    P_free(1:N) = 0;
                    Service(1:N) = 0;
                    Service2(1:N) = 0;
                        fprintf('=============================================\n');
                        fprintf('# Nodes: %d senders\n', N);
                        fprintf('CW: %d time slots\n', CW(1));
                        fprintf('Cycle: %d ms, duty-cycle ratio: %g\n', T, rho);
                        fprintf('Traffic: %g p/s\n', traffic_rate_bmac);
                        fprintf('=============================================\n');

                    %     for traffic_rate_bmac = [0.2 0.5 1 2] % X = X packets/s
                          for n = 1:N
                            out = model(sim_time, n, Ack_timeout_bamc, T, rho, CW(1), traffic_rate_bmac);
                            P_free(n) = out(1);
                            Pe(n) = out(2);
                            Service2(n) = out(3);
                            Service(n) = out(4);
                    %         fprintf('beta %i: %g\n', n, nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n));
                    %         SS = SS + nchoosek(N,n)*(1-P0(n))^n*P0(n)^(N-n)*S(n);
                          end

                        P0 = 0.01;

                        xguess=[P0]';
                        options = optimoptions('fsolve','Display','off'); %'iter'

                        xvect = fsolve('servtime_solver', xguess, options);

                        P0 = xvect(1);
                        fprintf('Expected P0: %g\n', P0);

                        ii = ii + 1;

                        SS = 0;
                        for i=1:N
                          SS = SS + nchoosek(N,i)*(1-P0)^i*P0^(N-i)*Pe(i);
                        end
                        thr(ii,jj) = 1-SS;
                        fprintf('Expected S: %g\n', thr(ii,jj));

    %                     SSS = 0;
    %                     for i=0:N-1
    %                       SSS = SSS + nchoosek(N-1,i)*(1-P0)^i*P0^(N-1-i)*Service2(i+1);
    %                     end
    % %                     ene(ii,jj) = (1-P0)*SSS/T;
    %                     ene(ii,jj) = SSS/T;

    %                     fprintf('Expected Servetime: %g\n', ((1-SS)*T/2 + SS*T));
    %                     
    %                     ene(ii,jj) = (1-P0)*((1-SS)*T/2 + SS*T)/T;
    % %                     ene(ii,jj) = (1-P0)*2*((1-SS)*T/2 + SS*T)/T;

                        ene(ii,jj) = (1-P0);
                        fprintf('Expected energy: %g\n', ene(ii,jj));

    %                     E_Service =(1-P0)/traffic_rate_bmac;
    %                     fprintf('E_Service: %g\n', E_Service*1e6);
    %                     ene(ii,jj) = (1-P0)*E_Service/T*1e6;
    %                     fprintf('Expected energy: %g\n', ene(ii,jj));






                    end
                end
            end
        end
%     end

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


