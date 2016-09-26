function y = tuning(x)

% clear all;
% clear all;
% clear all;


    global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
    global NN_B NN_W;
    global TR_B TR_W;
    global T_B T_W;
    global q0b q0w;
    global CASE;
    global TUNE_FACTOR TUNE_FACTOR2;
    global SS_W SS_B;
    
    global showed;
    
    


    %defaults:
    L_W = 150; %Packet size in DATA us
    Ackto_W = 10;

    cw_W = 16;
    CW_W(1:2) = [cw_W 16];%minimal CW size

    N_W = 40; % N_W/2 senders
    TR_W = 5;

    L_B = 1020;%3000;%us
    Ackto_B = 300;%256; %us

    cw_B = 70;
    CW_B(1:2) = [cw_B 70];

    N_B = 40; % N_W/2 senders
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
    
%     ii = 1;
%     jj = 1;

%     result_ThrB2{1,jj} = 'Thr_B';
%     result_ThrW2{1,jj} = 'Thr_W';
%     result_EneB2{1,jj} = 'Ene_B';
%     result_EneW2{1,jj} = 'Ene_W';



    NN_W = N_W/2;
    NN_B = N_B/2;

    NN_W_g = N_W/2; %for saved file generating
    NN_B_g = N_W/2; %for saved file generating


    T_W = round(x(1))*1e3; % ms
    Rho_W = round(x(2)*1e2)/1e2;
    T_B = round(x(3))*1e3; % ms
    Rho_B = round(x(4)*1e2)/1e2;

%     T_W = 100*1e3; % ms
%     T_B = 600*1e3; % ms
%     Rho_W = x(1);
%     Rho_B = x(2);


% T_W = 100*1e3;% ms
% Rho_W = round(0.07*1e2)/1e2;
% T_B = 401*1e3;% ms
% Rho_B = round(0.14*1e2)/1e2;


    fprintf('Trying T_W: %g, Rho_W: %g, T_B: %g, Rho_B: %g\n', T_W, Rho_W, T_B, Rho_B);
    
    outname=sprintf('result-(%g)-(%d-%g-%d)-(%d-%g-%d)-(%d-%d).mat', TR_B, ...
                            T_B, Rho_B, CW_B(1), ...
                            T_W, Rho_W, CW_W(1), TUNE_FACTOR, TUNE_FACTOR2);

%     load(outname);
    
    if ~exist(outname, 'file') %

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
% 
        save(outname, 'Sr_B', 'Pe_B', 'Sr_W', 'En_W', 'P_Wa', 'P_Wd');


    else
        fprintf('found the file\n');
        load(outname);
    end  

% return;    
    
%     continue;
    
    
%     jj = jj + 1;
%     result_ThrB2{1,jj} = sprintf('TR_B=%g mod', TR_B);
%     result_ThrW2{1,jj} = sprintf('TR_B=%g mod', TR_B);
%     result_EneB2{1,jj} = sprintf('TR_B=%g mod', TR_B);
%     result_EneW2{1,jj} = sprintf('TR_B=%g mod', TR_B);

%     ii = 1;
    
    if ~showed
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
        showed = 1;
    end

    [SS_B, SS_W, Q0_B, en_W] = compute(NN_B, TR_B, T_B, NN_W, TR_W, T_W);

%     ii = ii + 1;

%     result_ThrB2{ii,jj} = 1-SS_B;
%     result_ThrW2{ii,jj} = SS_W;
%     result_EneB2{ii,jj} = 1-Q0_B;
%     result_EneW2{ii,jj} = en_W;


%     result_ThrB2{ii,1} = sprintf('TR_W=%g', TR_W);
%     result_ThrW2{ii,1} = sprintf('TR_W=%g', TR_W);
%     result_EneB2{ii,1} = sprintf('TR_W=%g', TR_W);
%     result_EneW2{ii,1} = sprintf('TR_W=%g', TR_W);
    
    fprintf('Result Thr_W: %g, Thr_B: %g, Ene_W: %g, Ene_B: %g\n', SS_W, 1-SS_B, en_W, 1-Q0_B);
%     y = 1-Q0_B + en_W;
    y = (NN_B*(1-Q0_B) + NN_B*Rho_B + N_W*en_W)/(N_W+N_B);
   
    fprintf('Result: %g\n', y);
    
    

% end