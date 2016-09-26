
% function out = model( sim, no_B, to_B, cy_B, ra_B, cw_B, tr_B  )
function out = model(mod, nodes_B, packet_B, timeout_B, cycle_B, ratio_B, conwin_B, traffic_B, ...
                            nodes_W, packet_W, timeout_W, cycle_W, ratio_W, conwin_W)
% clear all;
% clear all;
% clear all;
% format short;


% nodes_B
% packet_B
% timeout_B
% cycle_B
% ratio_B
% conwin_B
% traffic_B
% 
% nodes_W
% packet_W
% timeout_W
% cycle_W
% ratio_W
% conwin_W




    globals;
    
    global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
    global TUNE_FACTOR TUNE_FACTOR2;
%     global PERCENTAGE;
%     global actualusPERCENTAGEage_d;
%     global energyusage_d;
%     global LL_Wd;
% for 
    no_B = nodes_B; % 10:5:30% 5:5:40 % %number of senders
    to_B = timeout_B; %timeout, us
    cy_B = cycle_B; %cycle length, us
    ra_B = ratio_B; %duty cycle ratio
    cw_B = conwin_B; %CW size
    tr_B = traffic_B;
    PERCENTAGE = min(max(tr_B/(1e6/cy_B),0.5),1);

%     PERCENTAGE = min(tr_B/(1e6/cy_B),1);

    no_W = nodes_W; %number of senders
    to_W = timeout_W; %timeout, us
    cy_W = cycle_W; %cycle length, us
    ra_W = ratio_W; %duty cycle ratio
    cw_W = conwin_W; %minimal CW size

    timeslot = 10; %us


    N_B = no_B;%senders
    tau_B = 2/3/(cw_B+10+1);%/4%1/40;
    Pa_B = ra_B;%1/10;
    Ds_B = packet_B;%3000;%us
    Df_B = Ds_B + to_B;%us
    % Df2 = Ds + timeout;%us
    Di_B = timeslot;%us
    T_B = cy_B; %1000000;
    Ta_B = T_B*Pa_B; %1000000;

    N_Wa = no_W;%senders
    Ds_W = 10;%150;%us ATIM
    Df_W = Ds_W + to_W;%us
    DIFS = 30;
    % Df2 = Ds + timeout;%us
    Di_W = timeslot;%us
    T_W = cy_W; %100 ms
    Ta_W = T_W*ra_W; %100 ms
    Pe_W = 1.0/Ta_W;
    W0_W = cw_W;
    m_W = 5;
    

    firstmactime_a(1:N_Wa+1,1:N_B+1) = 0;
    firstmactime_d(1:N_Wa+1,1:N_B+1) = 0;

    gamma_a(1:N_Wa+1,1:N_B+1) = 0;
    gamma_d(1:N_Wa+1,1:N_B+1) = 0;

    Sr_Ba(1:N_Wa+1,1:N_B+1) = 0;
    Sr_Bd(1:N_Wa+1,1:N_B+1) = 0;
%     Sr_B(1:N_Wa+1,1:N_B+1) = 0;

    actualusage_a(1:N_Wa+1,1:N_B+1) = 0;
    actualusage_d(1:N_Wa+1,1:N_B+1) = 0;

%     P_Wa(1:N_Wa+1,1:N_B+1) = 0;
%     P_Wd(1:N_Wa+1,1:N_B+1) = 0;

    Pe_Ba(1:N_Wa+1,1:N_B+1) = 0;
    Pe_Bd(1:N_Wa+1,1:N_B+1) = 0;
%     Pe_B(1:N_Wa+1,1:N_B+1) = 0;

%     Sr_W(1:N_Wa+1,1:N_B+1) = 0;


%     disp(TUNE_FACTOR*(max(N_B,15)/15-1) - TUNE_FACTOR2*alessthanb(N_B,15));
%     out = 1;
%     return;

    if DEBUG_HIGH; zw_printf('========= ATIM WINDOW =========\n'); end;


    for X_W = N_Wa:-1:0 %20%[20] % 21:40%20%

        % clear;
        % clear;
        % clear;
        % format long;




        % tau = 1/40;
        % 
        % P_awake = 1/10;
        % 
        % N = 40;
        % 
        % Ds = 3000;%us
        % Df1 = 300;%us
        % Df2 = 300;%us
        % Dns = 30;%us
        % T = 200000;


        Pi = 0.01;
        Pe_B_temp = 0.01;
        tau_W = 0.01;

        xguess=[Pi Pe_B_temp tau_W]';
        options = optimoptions('fsolve','Display','off'); %%'iter'
        xvect = fsolve('main_solver', xguess, options);
        % disp(xvect);

        Pi = xvect(1);
        Pe_B_temp = xvect(2);
        tau_W = xvect(3);

        Pe_Ba(X_W+1,N_B+1) = Pe_B_temp;

        if DEBUG_ATIM; zw_printf('\nPi(%d,%d): %g\n', N_B, X_W, Pi); end
        if DEBUG_ATIM; zw_printf('Pe_Ba(%d,%d): %g\n', N_B, X_W, Pe_B_temp); end
        if DEBUG_ATIM; zw_printf('tau_W(%d,%d): %g\n', N_B, X_W, tau_W); end

        out(1) = Pi;
        out(2) = Pe_B_temp;
        out(3) = tau_W;


    %     return;

        % % following is the new idea for LPL average MAC service time
        % %============================================================================================
        Sr_Ba(X_W+1,N_B+1) = (1-Pe_B_temp)*T_B*PERCENTAGE + Pe_B_temp*T_B;
    %     out(4) = Sr_Ba;
        if DEBUG_ATIM; zw_printf('Sr_Ba(%d,%d): %g\n', N_B, X_W, Sr_Ba(X_W+1,N_B+1)); end
        % %============================================================================================


        % % following is the new idea to get the number of successful nodes in ATIM/DATA window
        % %============================================================================================

    %     Ds_W = 10;


%         Pe_W = 0;

        temp_B = Pe_B_temp+(1-Pe_B_temp)*(1-tau_B);
        temp_W = Pe_W+(1-Pe_W)*(1-tau_W)^X_W;
        if DEBUG_ATIM; zw_printf('\ntemp_B(%d,%d): %g\n', N_B, X_W, temp_B); end
        if DEBUG_ATIM; zw_printf('temp_W(%d,%d): %g\n', N_B, X_W, temp_W); end


        %Channel model
        Pib0 = N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1)*Pa_B*temp_W;
        Pib1 = N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1)*(1-Pa_B)*temp_W;
        Pib2 = (1-temp_B^N_B-N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1))*temp_W;
        Pib3 = (1-Pe_W)*X_W*tau_W*(1-tau_W)^(X_W-1)*temp_B^N_B;
        Pib4 = ((1-Pe_W)*(1-(1-tau_W)^(X_W)-X_W*tau_W*(1-tau_W)^(X_W-1)))*temp_B^N_B;
        Pib5 = (1-temp_B^N_B)*(1-Pe_W)*(1-(1-tau_W)^X_W);
        if DEBUG_ATIM; zw_printf('\nPib0(%d,%d): %g\n', N_B, X_W, Pib0); end
        if DEBUG_ATIM; zw_printf('Pib1(%d,%d): %g\n', N_B, X_W, Pib1); end
        if DEBUG_ATIM; zw_printf('Pib2(%d,%d): %g\n', N_B, X_W, Pib2); end
        if DEBUG_ATIM; zw_printf('Pib3(%d,%d): %g\n', N_B, X_W, Pib3); end
        if DEBUG_ATIM; zw_printf('Pib4(%d,%d): %g\n', N_B, X_W, Pib4); end
        if DEBUG_ATIM; zw_printf('Pib5(%d,%d): %g\n', N_B, X_W, Pib5); end

        PiBs = Pib0; %A LPL transmission (successful)
        PiBc = Pib1 + Pib2 + Pib5; %A LPL transmission (collided (inter and intra))
        PiWs = Pib3; %successful PSM transmission
        PiWc = Pib4; %intra-collided PSM transmission
        Pii = 1-PiBs-PiBc-PiWc-PiWs; %no transmission
        PiWns = 1 - PiWs; %Not a successful PSM transmission

        if DEBUG_ATIM; zw_printf('\nPiBs(%d,%d)=Pib0: %g\n', N_B, X_W, PiBs); end
        if DEBUG_ATIM; zw_printf('PiBc(%d,%d)=Pib1+Pib2+Pib5: %g\n', N_B, X_W, PiBc); end
        if DEBUG_ATIM; zw_printf('PiWs(%d,%d)=Pib3: %g\n', N_B, X_W, PiWs); end
        if DEBUG_ATIM; zw_printf('PiWc(%d,%d)=Pib4: %g\n', N_B, X_W, PiWc); end
        if DEBUG_ATIM; zw_printf('Pii(%d,%d)=1-PiB-PiWc-PiWs: %g\n', N_B, X_W, Pii); end
        if DEBUG_ATIM; zw_printf('PiWns(%d,%d)=1-PiWs: %g\n', N_B, X_W, PiWns); end


    %     ALf_W = (Pii*timeslot + PiBs/(1-Pii)*Ds_B + PiBc/(1-Pii)*Df_B + PiWc/(1-Pii)*(Ds_W+to_W+DIFS))/PiWns;%;%
        ALf_W = (Pii*timeslot + PiBs*Ds_B + PiBc*Df_B + PiWc*(Ds_W+to_W+DIFS))/PiWns;%;%
        if DEBUG_ATIM; zw_printf('\nALf_W(%d,%d): %g\n', N_B, X_W, ALf_W); end

    %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs+N_B/10) + Ds_W;
        if N_B~=0
            Sr_first_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs) + Ds_W + ...
               TUNE_FACTOR*(1-alessthanb(N_B,10))*(N_B/10-1) - TUNE_FACTOR2*alessthanb(N_B,10);% TUNE_FACTOR*(max(N_B,10)/10-1);% 
        else
            Sr_first_W = ALf_W * (PiWns/PiWs) + Ds_W;%TUNE_FACTOR*(N_B/10-1);%
        end
            %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs) + Ds_W + DIFS + 1/tau_W*timeslot + 400*(N_B/10-1);
    %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs+N_B/10) + Ds_W + DIFS + 1/tau_W*timeslot;
        if DEBUG_ATIM; zw_printf('Sr_first_Wa(%d,%d): %g\n', N_B, X_W, Sr_first_W); end
%         fprintf('Sr_first_Wa(%d,%d): %g\n', N_B, X_W, Sr_first_W);
    %     firstmactime(X_W,N_B/5) = Sr_W;
    %     if Sr_W == Inf
    %         Sr_W = 1e6;
    %     end

        firstmactime_a(X_W+1,N_B+1) = Sr_first_W;



        % fprintf('(1-Pi)*Ds_B/2: %g\n', (1-Pi)*Ds_B/2);
        % fprintf('PiWns/PiWs: %g\n', PiWns/PiWs);

        % sumsum = 0;
        % for n=1:15
        %     for x=1:n
        %         
        %         
        %         for y=1:n-x
        %           
        %             sumsum = sumsum + nchoosek(n,x)*Pii^x*nchoosek(n-x,y)*PiB^y*PiWc^(n-x-y)*PiWs*(x+y*(Ds_B)+(n-x-y+1)*Ds_W);
        %             
        %         end
        %         
        %         
        %     end
        %     
        % end
        % 
        % 
        % Sr_W = (1-Pi)*Ds_B/2 + sumsum;
        % fprintf('\nSr_W(%d,%d): %g\n', N_B, X_W, Sr_W);

        % %============================================================================================


    end



    Actnode_stat(1:T_W/timeslot)=0; % 10us per unit
    index = 1;

    for i=N_Wa:-1:0

        if i==0 %firstmactime(i+1,N_B+1) == Inf
            firstmactime_a(i+1,N_B+1) = T_W*ra_W - index*timeslot;
        end

        Actnode_stat(index:index+round(firstmactime_a(i+1,N_B+1)/timeslot) )= i;
%         i
%         firstmactime_a(i+1,N_B+1)
%         index+round(firstmactime_a(i+1,N_B+1)/timeslot)
        index = index + round(firstmactime_a(i+1,N_B+1)/timeslot);
        actualusage_a(i+1,N_B+1) = actualusage_a(i+1,N_B+1)+round(firstmactime_a(i+1,N_B+1)/timeslot);

        if index > T_W*ra_W/timeslot%T_W/timeslot %

            if i>0 %firstmactime(i+2,N_B+1) ~= Inf
                temp = index - T_W*ra_W/timeslot;
                actualusage_a(i+1,N_B+1) = actualusage_a(i+1,N_B+1)- temp+1;
            end

            break;
        end

    end
    
    

% T_W/timeslot
% actualusage_a

    gamma_a=actualusage_a/(T_W/timeslot);



%     pause;


    N_Wd = N_Wa-i;
    P_Wa(N_Wa+1,N_B+1) = N_Wd/N_Wa;
    
%     if N_Wd==0 && N_Wa==1
%         N_Wd = 1;
%         P_Wa(N_Wa+1,N_B+1) = N_Wd/N_Wa;
%     end

    Ds_W = packet_W;%us
    Df_W = Ds_W + to_W;%us

    if DEBUG_HIGH; zw_printf('\n\n========= DATA WINDOW =========\n'); end;
    Pe_W = 1/T_W;


    for X_W = N_Wd:-1:0 %20%[20] % 21:40%20%

        % clear;
        % clear;
        % clear;
        % format long;




        % tau = 1/40;
        % 
        % P_awake = 1/10;
        % 
        % N = 40;
        % 
        % Ds = 3000;%us
        % Df1 = 300;%us
        % Df2 = 300;%us
        % Dns = 30;%us
        % T = 200000;


        Pi = 0.01;
        Pe_B_temp = 0.01;
        tau_W = 0.1;

        xguess=[Pi Pe_B_temp tau_W]';
        options = optimoptions('fsolve','Display','off'); %%'iter'
        xvect = fsolve('main_solver', xguess, options);
        % disp(xvect);

        Pi = xvect(1);
        Pe_B_temp = xvect(2);
        tau_W = xvect(3);

        if DEBUG_DATA; zw_printf('\nPi(%d,%d): %g\n', N_B, X_W, Pi); end
        if DEBUG_DATA; zw_printf('Pe_Bd(%d,%d): %g\n', N_B, X_W, Pe_B_temp); end
        if DEBUG_DATA; zw_printf('tau_W(%d,%d): %g\n', N_B, X_W, tau_W); end

        out(1) = Pi;
        out(2) = Pe_B_temp;
        out(3) = tau_W;

        Pe_Bd(X_W+1,N_B+1) = Pe_B_temp;


        % return;

        % % following is the new idea for LPL average MAC service time
        % %============================================================================================

        Sr_Bd(X_W+1,N_B+1) = (1-Pe_B_temp)*T_B*PERCENTAGE + Pe_B_temp*T_B;
    %     out(4) = Sr_Bd;
        if DEBUG_DATA; zw_printf('Sr_Bd(%d,%d): %g\n', N_B, X_W, Sr_Bd(X_W+1,N_B+1)); end
        % %============================================================================================


        % % following is the new idea to get the number of successful nodes in ATIM/DATA window
        % %============================================================================================

    %     Ds_W = 10;


%         Pe_W = 0;

        temp_B = Pe_B_temp+(1-Pe_B_temp)*(1-tau_B);
        temp_W = Pe_W+(1-Pe_W)*(1-tau_W)^X_W;
        if DEBUG_DATA; zw_printf('\ntemp_B(%d,%d): %g\n', N_B, X_W, temp_B); end
        if DEBUG_DATA; zw_printf('temp_W(%d,%d): %g\n', N_B, X_W, temp_W); end


        %Channel model
        Pib0 = N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1)*Pa_B*temp_W;
        Pib1 = N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1)*(1-Pa_B)*temp_W;
        Pib2 = (1-temp_B^N_B-N_B*tau_B*(1-Pe_B_temp)*temp_B^(N_B-1))*temp_W;
        Pib3 = (1-Pe_W)*X_W*tau_W*(1-tau_W)^(X_W-1)*temp_B^N_B;
        Pib4 = ((1-Pe_W)*(1-(1-tau_W)^(X_W)-X_W*tau_W*(1-tau_W)^(X_W-1)))*temp_B^N_B;
        Pib5 = (1-temp_B^N_B)*(1-Pe_W)*(1-(1-tau_W)^X_W);
        if DEBUG_DATA; zw_printf('\nPib0(%d,%d): %g\n', N_B, X_W, Pib0); end
        if DEBUG_DATA; zw_printf('Pib1(%d,%d): %g\n', N_B, X_W, Pib1); end
        if DEBUG_DATA; zw_printf('Pib2(%d,%d): %g\n', N_B, X_W, Pib2); end
        if DEBUG_DATA; zw_printf('Pib3(%d,%d): %g\n', N_B, X_W, Pib3); end
        if DEBUG_DATA; zw_printf('Pib4(%d,%d): %g\n', N_B, X_W, Pib4); end
        if DEBUG_DATA; zw_printf('Pib5(%d,%d): %g\n', N_B, X_W, Pib5); end

        PiBs = Pib0; %A LPL transmission (successful)
        PiBc = Pib1 + Pib2 + Pib5; %A LPL transmission (collided (inter and intra))
        PiWs = Pib3; %successful PSM transmission
        PiWc = Pib4; %intra-collided PSM transmission
        Pii = 1-PiBs-PiBc-PiWc-PiWs; %no transmission
        PiWns = 1 - PiWs; %Not a successful PSM transmission

        if DEBUG_DATA; zw_printf('\nPiBs(%d,%d)=Pib0: %g\n', N_B, X_W, PiBs); end
        if DEBUG_DATA; zw_printf('PiBc(%d,%d)=Pib1+Pib2+Pib5: %g\n', N_B, X_W, PiBc); end
        if DEBUG_DATA; zw_printf('PiWs(%d,%d)=Pib3: %g\n', N_B, X_W, PiWs); end
        if DEBUG_DATA; zw_printf('PiWc(%d,%d)=Pib4: %g\n', N_B, X_W, PiWc); end
        if DEBUG_DATA; zw_printf('Pii(%d,%d)=1-PiB-PiWc-PiWs: %g\n', N_B, X_W, Pii); end
        if DEBUG_DATA; zw_printf('PiWns(%d,%d)=1-PiWs: %g\n', N_B, X_W, PiWns); end


    %     ALf_W = (Pii*timeslot + PiBs/(1-Pii)*Ds_B + PiBc/(1-Pii)*Df_B + PiWc/(1-Pii)*(Ds_W+to_W+DIFS))/PiWns;%;%
        ALf_W = (Pii*timeslot + PiBs*Ds_B + PiBc*Df_B + PiWc*(Ds_W+to_W+DIFS))/PiWns;%;%
        if DEBUG_DATA; zw_printf('\nALf_W(%d,%d): %g\n', N_B, X_W, ALf_W); end

    %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs+N_B/10) + Ds_W;
        if N_B~=0
            Sr_first_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs) + Ds_W  + ...
                TUNE_FACTOR*(1-alessthanb(N_B,10))*(N_B/10-1) - TUNE_FACTOR2*alessthanb(N_B,10);% TUNE_FACTOR*(max(N_B,10)/10-1);%
        else
            Sr_first_W = ALf_W * (PiWns/PiWs) + Ds_W;% TUNE_FACTOR*(max(N_B,10)/10-1);%
        end

        %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs) + Ds_W + DIFS + 1/tau_W*timeslot + 400*(N_B/10-1);
    %     Sr_W = (1-Pi)/(PiBs+PiBc)*(PiBc*Df_B/2 + PiBs*Ds_B/2) + ALf_W * (PiWns/PiWs+N_B/10) + Ds_W + DIFS + 1/tau_W*timeslot;
        if DEBUG_DATA; zw_printf('Sr_first_Wd(%d,%d): %g\n', N_B, X_W, Sr_first_W); end
    %     firstmactime(X_W,N_B/5) = Sr_W;
    %     if Sr_W == Inf
    %         Sr_W = 1e6;
    %     end
        firstmactime_d(X_W+1,N_B+1) = Sr_first_W;


        % fprintf('(1-Pi)*Ds_B/2: %g\n', (1-Pi)*Ds_B/2);
        % fprintf('PiWns/PiWs: %g\n', PiWns/PiWs);

        % sumsum = 0;
        % for n=1:15
        %     for x=1:n
        %         
        %         
        %         for y=1:n-x
        %           
        %             sumsum = sumsum + nchoosek(n,x)*Pii^x*nchoosek(n-x,y)*PiB^y*PiWc^(n-x-y)*PiWs*(x+y*(Ds_B)+(n-x-y+1)*Ds_W);
        %             
        %         end
        %         
        %         
        %     end
        %     
        % end
        % 
        % 
        % Sr_W = (1-Pi)*Ds_B/2 + sumsum;
        % fprintf('\nSr_W(%d,%d): %g\n', N_B, X_W, Sr_W);

        % %============================================================================================



    end

    index = round(T_W*ra_W/timeslot+1);
%     fprintf('%d\n', index);
%     energyusage_d = 0;
    j = 1;
    energyusage_d(j) = Inf;
    for i=N_Wd:-1:0

%         Actnode_stat(index:index+round(Sr_W) )= i;
%         round(Sr_W )
%         round(firstmactime(i+1,1)/timeslot)
        if i==0 %firstmactime(i+1,N_B+1) == Inf
            firstmactime_d(i+1,N_B+1) = T_W - index*timeslot;
        end
%         try
%             fprintf('%g, ', round(firstmactime_d(i+1,N_B+1)/timeslot));
%             fprintf('%d\n', index);
            Actnode_stat(index:index+round((firstmactime_d(i+1,N_B+1)/timeslot)) )= i;
%         catch ME
%             ME.identifier
%             error('')
%         end
        index = index + round(firstmactime_d(i+1,N_B+1)/timeslot);
        actualusage_d(i+1,N_B+1) = actualusage_d(i+1,N_B+1) + round(firstmactime_d(i+1,N_B+1)/timeslot);
        if i~=0
            energyusage_d(j) = round(sum(firstmactime_d(N_Wd+2-j:N_Wd+1,N_B+1)));
    %         N_Wd+2-j
    %         firstmactime_d(N_Wd+2-j:N_Wd+1,N_B+1)
            j = j + 1;
        end
        if index > T_W/timeslot%T_W/timeslot %

            if i>0 %firstmactime(i+2,N_B+1) ~= Inf
                temp = index - T_W/timeslot;
                actualusage_d(i+1,N_B+1) = actualusage_d(i+1,N_B+1)- temp+1;
            end

            break;
        end

    end    
    
%     disp(energyusage_d);
% 
    P_Wd(N_Wa+1,N_B+1) = (N_Wd-i)/N_Wd;
%     N_Wa
    if N_Wa~=0
        LL_Wd(N_Wa+1,N_B+1) = Ta_W + sum(energyusage_d)/N_Wd;
%         disp(LL_Wd(N_Wa+1,N_B+1));
%         disp(Ta_W);
    
%     Ta_W/T_W
        if P_Wa(N_Wa+1,N_B+1)~=1
            En_W(N_Wa+1,N_B+1) = Ta_W/T_W + ((LL_Wd(N_Wa+1,N_B+1)-Ta_W)*P_Wa(N_Wa+1,N_B+1)/T_W/(1-P_Wa(N_Wa+1,N_B+1)))*-log(P_Wa(N_Wa+1,N_B+1));
%             disp(P_Wa(N_Wa+1,N_B+1)/(1-P_Wa(N_Wa+1,N_B+1))*-log(P_Wa(N_Wa+1,N_B+1)));
        else
            En_W(N_Wa+1,N_B+1) = LL_Wd(N_Wa+1,N_B+1)/T_W;
        end
    end
    
%     disp(En_W(N_Wa+1,N_B+1));
%     (LL_Wd(N_Wa+1,N_B+1)*P_Wa(N_Wa+1,N_B+1)/T_W/(1-P_Wa(N_Wa+1,N_B+1)))
%     P_Wa(N_Wa+1,N_B+1)
%     (LL_Wd(N_Wa+1,N_B+1)*P_Wa(N_Wa+1,N_B+1)/T_W/(1-P_Wa(N_Wa+1,N_B+1)))*log(P_Wa(N_Wa+1,N_B+1))
%     -log(P_Wa(N_Wa+1,N_B+1))
    
%     ssss=0;
%     q=0.7;
%     for i=1:1000
%         ssss = ssss + q^i/i;
%         
%     end
%     ssss
%     log(1-q)

    gamma_d=actualusage_d/(T_W/timeslot);

%     asfig = figure;
%     plot(Actnode_stat, 'DisplayName','Actnode stat');
%     xlabel('Time');ylabel('Average # Active nodes A(t)');
%     title(sprintf('N_W=%d, W_W=%d, T_W=%d us, L_W=%d, N_B=%d, W_B=%d, T_B=%d us, L_B=%d', N_Wa, cw_W, T_W/1e3, Ds_W, N_B, cw_B, T_B/1e3, Ds_B));

    if DEBUG_HIGH; zw_printf('\n\n========= OVERALL =========\n'); end;

    Sr_B(N_Wa+1,N_B+1) = gamma_a(:,N_B+1)'*Sr_Ba(:,N_B+1) + gamma_d(:,N_B+1)'*Sr_Bd(:,N_B+1);
%     fprintf('Sr_Ba(%d,%d): %g\n', N_B, N_Wa, Sr_Ba(N_Wa+1,N_B+1));
%     fprintf('Sr_Bd(%d,%d): %g\n', N_B, N_Wa, Sr_Bd(N_Wa+1,N_B+1));
    if DEBUG_HIGH; zw_printf('Sr_B(%d,%d): %g\n', N_B, N_Wa, Sr_B(N_Wa+1,N_B+1)); end;

    Pe_B(N_Wa+1,N_B+1) = gamma_a(:,N_B+1)'*Pe_Ba(:,N_B+1) + gamma_d(:,N_B+1)'*Pe_Bd(:,N_B+1);
%     fprintf('Pe_Ba(%d,%d): %g\n', N_B, N_Wa, Pe_Ba(N_Wa+1,N_B+1));
%     fprintf('Pe_Bd(%d,%d): %g\n', N_B, N_Wa, Pe_Bd(N_Wa+1,N_B+1));
    if DEBUG_HIGH; zw_printf('Pe_B(%d,%d): %g\n', N_B, N_Wa, Pe_B(N_Wa+1,N_B+1)); end;

%     for i = 0:100
%         Sr_W(N_Wa+1,N_B+1) = Sr_W(N_Wa+1,N_B+1) + (i+1)*T_W*(1-P_Wa(N_Wa+1,N_B+1)*P_Wd(N_Wa+1,N_B+1))^i;
%     end
    Sr_W(N_Wa+1,N_B+1) = T_W/(P_Wa(N_Wa+1,N_B+1)*P_Wd(N_Wa+1,N_B+1));
    if DEBUG_HIGH; zw_printf('Sr_W(%d,%d): %g\n', N_B, N_Wa, Sr_W(N_Wa+1,N_B+1)); end;
    if DEBUG_HIGH; zw_printf('En_W(%d,%d): %g\n', N_B, N_Wa, En_W(N_Wa+1,N_B+1)); end;
    if DEBUG_HIGH; zw_printf('P_Wa(%d,%d): %g\n', N_B, N_Wa, P_Wa(N_Wa+1,N_B+1)); end;


% end
end





