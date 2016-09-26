function [] = wifi_psm( TT, n, payload, tr, atim, data )

    rng('default');
    rng(0);

    %nbFailure(1:N)=0;
    %nbCCA(1:N)=0;

    STATE_IDLE=0; % for waiting for a packet
    STATE_SOFT=1; % for software delay
    STATE_INIT=2;
    STATE_DIFS=3;
    STATE_CSMA=4;
    STATE_TX=5;
    STATE_WAIT=6; % for PSM only, means waiting for ATIM/DATA/BEACON window to come.

    cnt=0;
    rd=0;

    PAYLOAD = payload;%bytes
    HEDEAR = 40;%bytes, mac header, tcp/udp header and ip header
    PREAMBLE = 4;%40 us
    aMinBE=4;
    aMaxBE=aMinBE+5;
    LDATA=PREAMBLE+ceil(HEDEAR*8/54/10)+ceil(PAYLOAD*8/54/10);%ceil(846*8/54/10)
    ATIM=30;%bytes
    LATIM=PREAMBLE+ceil(ATIM*8/54/10);
    LACK=0;
    LACKTIMEOUT=0;
    DIFS=1;


    % disp(LDATA);
    % disp(LATIM);
    % return;


    T=TT;
    queue_size = 100;
    SOFT_DELAY = 1; %assume the software delay is very small
    SATURATION_ENABLED = 0; %0: unsaturated
    TICK_INTERVAL = 50000; %show interval
    DYNAMIC_PLOT = 1; %draw the queue size gragh dynamically
    UPDATE_INTERVAL = 10000; %drawing interval
    SHOW_INTERVAL = 500000; %drawing interval
    traffic_rate = tr; %packets/s

    Trials=n;%number of nodes

    %window length
    BEACON_WINDOW = 1; % *10 us
    ATIM_WINDOW = atim; % *10 us
    DATA_WINDOW = data; % *10 us
    TOTAL_WINDOW = BEACON_WINDOW + ATIM_WINDOW + DATA_WINDOW;

    fprintf('# Nodes: %d\n', Trials);
    fprintf('Simulation time: %g s\n', T/100000);
    fprintf('CW: %d, max: %d, limit: N/A\n', 2^aMinBE, 2^aMaxBE);
    fprintf('====ATIM window related:====\n');
    fprintf('packet size: %d, window size: %d\n', LATIM, ATIM_WINDOW);
    fprintf('====DATA window related:====\n');
    fprintf('packet size: %d, DATA window: %d\n', LDATA, DATA_WINDOW);

    fprintf('lambda: %f\n', traffic_rate);
    fprintf('mu: < %f\n', 1000/((ATIM_WINDOW+DATA_WINDOW)/100));



    %time window
    TIMESTATE_BEACON = 0; %beacon window
    TIMESTATE_ATIM = 1; %ATIM window
    TIMESTATE_DATA = 2; %DATA window


    %debug level
    %===========================
    debug_backoff = 0;
    debug_tx = 0;
    debug_queue = 0;
    debug_psm = 0;
    debug_collision = 0;

    %===========================



    % zw_time(1:6)=0;
    % zw_time=clock();
    % outname=['sim-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
    %     num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
    %     '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];
    % diary(outname);
    % diary on;

    for N=Trials:10:Trials
        %disp(N);

        averageCCA(1:Trials)=0;
        averageBusy1(1:Trials)=0;
        averageBusy2(1:Trials)=0;
        averageTrans(1:Trials)=0;
        totalTrans(1:Trials)=0;
        averageSucTrans(1:Trials)=0;
        collisions(1:Trials)=0;
        delay(1:Trials)=0;
        CW(1:Trials)=aMinBE;% aMinBE < BE < aMaxBE
        busyFor(1:Trials)=0; %for transmission delay
        tempbusyFor(1:Trials)=0; %for transmission delay
        nbTransmission(1:Trials)=0;
        nbDataTransmission(1:Trials)=0;%Successful tx
        nbCollision(1:Trials)=0;
        successTX(1:Trials)=2;%0 fail, 1 success, 2 no data
        DIFScnt(1:Trials)=0; %for DIFS delay
        oneshot(1:Trials)=0; %for oneshot of transmission, 1 means the first try, 0 means the first try fails.
        if SATURATION_ENABLED == 1
            nodeState(1:Trials)=STATE_WAIT; %for FSM 
        else
            nodeState(1:Trials)=STATE_WAIT; %for FSM 
        end
        trfail(1:Trials)=0; %for tr fails.
        softdelay(1:Trials)=SOFT_DELAY;
        interarrival(1:Trials)=-1;
        
        TXtime(1:Trials)=0;
        Awaketime(1:Trials)=0;
        LastBeacontime(1:Trials) = 0;
        StartTimeofCurrentCycle = 0;

        stateintime = TIMESTATE_BEACON;



        mac_queue(1:Trials)=queue_size; %current queue residual
        queue_empty_count(1:Trials)=0;
        queue_size_count(1:Trials)=0;

        if DYNAMIC_PLOT == 1

            t=[0];
            m=[];
            for i=1:N
                m=[m;0];
            end

            p = plot(t,m,...
               'EraseMode','background','MarkerSize',5);hold on;
            axis([1 SHOW_INTERVAL/UPDATE_INTERVAL -5 queue_size+5]);
            grid off;

        end    


        tstart = clock;

        global_sync = 0;

        tt=tic;
        for i=1:T
            if mod(i,TICK_INTERVAL)==0
                diary off;
                fprintf('%gs\n', i/100000);
                toc(tt);
                diary on;
                tt=tic;
            end


            if DYNAMIC_PLOT == 1

                %draw figure to show the queue usage
                if mod(i,UPDATE_INTERVAL)==0
                    t=[t i/UPDATE_INTERVAL];                   %Matrix 1*(i+1)
    %                     disp(size(t,1));
    %                     disp(size(t,2));
                    if size(t,2)==SHOW_INTERVAL/UPDATE_INTERVAL+2
                        t=t(1:SHOW_INTERVAL/UPDATE_INTERVAL+1);
    %                     disp(size(t,1));
    %                     disp(size(t,2));
    %                     return;
                    end
                    n=[];
                end
            end

    %         fprintf('\n');

            for j=1:N

                if DYNAMIC_PLOT == 1
                    if mod(i,UPDATE_INTERVAL)==0
                        n=[n;queue_size-mac_queue(j)];
                    end
                end

                %Poisson traffic generator
                %===============================================================================
                if SATURATION_ENABLED == 0
                    
                    queue_size_count(j) = queue_size_count(j) + queue_size-mac_queue(j);

                    if interarrival(j)==-1
                        while(1)
                            interarrival(j) = ceil((-1/traffic_rate)*log(rand)*1e5);%ceil(exprnd(1/traffic_rate)*1e5);
%                             fprintf('[%d us @ %d 1. interarrival: %d us\n]', (i-1)*10, j, interarrival(j)*10);
                            if interarrival(j)==0
                                %put data into mac_queue
                                if mac_queue(j)==0
                                    if debug_queue; fprintf('[%d us @ %d] !!!!!!no space in the mac_queue, drop it!!!!!\n', (i-1)*10, j); end
                                else
                                    mac_queue(j) = mac_queue(j) - 1;
    %                                 fprintf('[%d us] put data into the mac_queue(%d)\n', (i-1)*10, j);
                                end
                            else
                                break;
                            end
                        end
                        interarrival(j) = interarrival(j) - 1;
                        if interarrival(j)==0 
        %                     fprintf('################# interarrival(%d): %d !!!!!!!!!!!!!!!!\n', j, interarrival(j)*10);
                            interarrival(j)=1;
                        end
    %                     fprintf('[%d us] interarrival(%d): %d\n', (i-1)*10, j, interarrival(j)*10);

                    elseif interarrival(j)==1
                        %put data into mac_queue
                        if mac_queue(j)==0
                            if debug_queue; fprintf('[%d us @ %d] no space in the mac_queue, drop it!!!!!\n', (i-1)*10, j); end
                        else
                            mac_queue(j) = mac_queue(j) - 1;
    %                         fprintf('[%d us] put data into the mac_queue(%d)\n', (i-1)*10, j);
                        end

                        while(1)
                            interarrival(j) = ceil((-1/traffic_rate)*log(rand)*1e5);%ceil(exprnd(1/traffic_rate)*1e5);
%                             fprintf('[%d us @ %d 2. interarrival: %d us\n]', (i-1)*10, j, interarrival(j)*10);
                            if interarrival(j)==0
                            %put data into mac_queue
                                if mac_queue(j)==0
                                    if debug_queue; fprintf('[%d us @ %d] no space in the mac_queue, drop it!!!!!\n', (i-1)*10, j); end
                                else
                                    mac_queue(j) = mac_queue(j) - 1;
    %                                 fprintf('[%d us] put data into the mac_queue(%d)\n', (i-1)*10, j);
                                end
                            else
                                break;
                            end
                        end

                        interarrival(j) = interarrival(j) - 1;
                        if interarrival(j)==0
        %                     fprintf('################# interarrival(%d): %d !!!!!!!!!!!!!!!!\n', j, interarrival(j)*10);
                            interarrival(j)=1;
                        end
    %                     fprintf('[%d us] interarrival(%d): %d\n', (i-1)*10, j, interarrival(j)*10);
                    else
                        interarrival(j) = interarrival(j) - 1;

                    end

                    if mac_queue(j)==queue_size%queue is empty

                        queue_empty_count(j) = queue_empty_count(j) + 1;

                    end



                end
                %===============================================================================


                if global_sync == BEACON_WINDOW % BEACON WINDOW is over
                    %prepare to send ATIM if j has data to send
                    StartTimeofCurrentCycle = i;
                    if SATURATION_ENABLED == 0
                        if successTX(j)==1 || successTX(j)==2 %send successfully or no data in last round, pick a new packet in the queue.
                            if mac_queue(j)~=queue_size
                                %fetch data from mac_queue
%                                 mac_queue(j) = mac_queue(j) + 1;
            %                     fprintf('[%d us] get a data in mac_queue(%d)\n', (i-1)*10, j);
                                LastBeacontime(j) = i;
                                nodeState(j)=STATE_INIT;
                                successTX(j)=0;
                                if debug_psm; fprintf('[%d us @ %d] beacon window ends, begin to send ATIM... mac_queue: %d\n', (i-1)*10, j, mac_queue(j)); end
                            else %no data to send
                                nodeState(j)=STATE_WAIT;
                                successTX(j)=2;
                                if debug_psm; fprintf('[%d us @ %d] beacon window ends, but no data to send...\n', (i-1)*10, j); end
                            end

                        else % send fail in last round, still use the old packet.
                            if debug_psm; fprintf('[%d us @ %d] beacon window ends, begin to send ATIM (old)... mac_queue: %d\n', (i-1)*10, j, mac_queue(j)); end
                            nodeState(j)=STATE_INIT;
                        end
                    else
                        nodeState(j)=STATE_INIT;
                        successTX(j)=0;
                        if debug_psm; fprintf('[%d us @ %d] beacon window ends, begin to send ATIM...\n', (i-1)*10, j); end
                    end
                    stateintime = TIMESTATE_ATIM;

                elseif global_sync == BEACON_WINDOW + ATIM_WINDOW % ATIM WINDOW is over
                    %check if a node has sent an ATIM successfully.
    %                 if debug_psm; fprintf('[%d us @ %d] ATIM window ends...\n', (i-1)*10, j); end
                    stateintime = TIMESTATE_DATA;

                    %if yes, then set the node state to STATE_INIT
                    if successTX(j)==1
                        if debug_psm; fprintf('[%d us @ %d] ATIM sends OK, tx data now...\n', (i-1)*10, j); end
                        nodeState(j)=STATE_INIT; %begin to transmit data
                        successTX(j)=0;
    %                     if debug_psm; fprintf('[%d us @ %d] delay: %d, busyFor: %d, tempbusyFor: %d, DIFScnt: %d, trfail: %d\n', ...
    %                         (i-1)*10, j, delay(j), busyFor(j), tempbusyFor(j), DIFScnt(j), trfail(j)); end


                    %otherwise, do nothing.
                    elseif successTX(j)==0
                        if debug_psm; fprintf('[%d us @ %d] ATIM sends not OK, retry in next round... awaketime: %f += %f\n', (i-1)*10, j, Awaketime(j), BEACON_WINDOW + ATIM_WINDOW); end
                        Awaketime(j) = Awaketime(j) + (BEACON_WINDOW + ATIM_WINDOW);
                        nodeState(j)=STATE_WAIT;
                        successTX(j)=0; 
                        delay(j)=0; busyFor(j)=0; tempbusyFor(j)=0; DIFScnt(j)=0; trfail(j)=0;
    %                     if debug_psm; fprintf('[%d us @ %d] delay: %d, busyFor: %d, tempbusyFor: %d, DIFScnt: %d, trfail: %d\n', ...
    %                         (i-1)*10, j, delay(j), busyFor(j), tempbusyFor(j), DIFScnt(j), trfail(j)); end
                    else
                        if debug_psm; fprintf('[%d us @ %d] no data to send...\n', (i-1)*10, j); end
                    end


                elseif global_sync == TOTAL_WINDOW % DATA WINDOW is over
                    if j==N
                        global_sync = 0;
                    end
                    stateintime = TIMESTATE_BEACON;
    %                 if debug_psm; fprintf('[%d us @ %d] DATA window ends...\n', (i-1)*10, j); end

                    if successTX(j)==1
%                         if debug_psm; fprintf('[%d us @ %d] DATA sends OK, send new ATIM+DATA in next round... awaketime: %f += %f\n', (i-1)*10, j, Awaketime(j), i - StartTimeofCurrentCycle); end
                        nbDataTransmission(j) = nbDataTransmission(j)+1;
                        TXtime(j) = TXtime(j) + (i - LastBeacontime(j) + BEACON_WINDOW);
                        mac_queue(j) = mac_queue(j) + 1;
    %                     successTX(j)=0;

                    %otherwise, do nothing.
                    elseif successTX(j)==0
                        if debug_psm; fprintf('[%d us @ %d] DATA sends not OK, retry in next round...', (i-1)*10, j); end
    %                     successTX(j)=0; 
                        if nodeState(j)~=STATE_WAIT
                            if debug_psm; fprintf('awaketime: %f += %f\n', Awaketime(j), TOTAL_WINDOW);end
                            Awaketime(j) = Awaketime(j) + (TOTAL_WINDOW);
                        else
                            if debug_psm; fprintf('\n');end
                        end
                        delay(j)=0; busyFor(j)=0; tempbusyFor(j)=0; DIFScnt(j)=0; trfail(j)=0;
                    else
                        if debug_psm; fprintf('[%d us @ %d] no data to send...\n', (i-1)*10, j); end
                    end
                    nodeState(j)=STATE_WAIT;
                end


%                 if nodeState(j)==STATE_IDLE
%                     %check if mac_queue has data
% 
%                     if mac_queue(j)~=queue_size
%                         %fetch data from mac_queue
%                         mac_queue(j) = mac_queue(j) + 1;
%     %                     fprintf('[%d us] get a data in mac_queue(%d)\n', (i-1)*10, j);
%                         nodeState(j)=STATE_SOFT;
%                     else
%                     end
%                 else
                if nodeState(j)==STATE_WAIT

                    %do nothing
    %                 fprintf('[%d us] waiting...\n', (i-1)*10);
                    continue;

                elseif nodeState(j)==STATE_INIT
                   DIFScnt(j)=0;
                    if trfail(j)==0
                        oneshot(j)=0;
                    else
                        oneshot(j)=0;
                    end
                    nodeState(j)=STATE_DIFS;
                    continue;
                end

                if nodeState(j)==STATE_SOFT
                    softdelay(j)=softdelay(j)-1;
                    if softdelay(j)==0
                        nodeState(j)=STATE_INIT;
                        softdelay(j)=SOFT_DELAY;
                    end
                end

                switch nodeState(j)
        %             case{STATE_INIT}
        %                 DIFScnt(j)=0;
        %                 oneshot(j)=1;
        %                 nodeState(j)=STATE_DIFS;
                    case{STATE_DIFS}
                        if sum(busyFor(1:N))==0 % channel idle
                            %disp('channel idle');
                            if debug_backoff; fprintf('[%d us @ %d] channel is idle, current state: STATE_DIFS.\n',(i-1)*10, j); end
                            DIFScnt(j)=DIFScnt(j)+1;
                        else % channle not idle
                            if debug_backoff; fprintf('[%d us @ %d] channel is busy, current state: STATE_DIFS.\n',(i-1)*10, j); end
                            oneshot(j)=0;
                            DIFScnt(j)=0;
                        end
                    case{STATE_CSMA}
                        if sum(busyFor(1:N))==0 % channel idle
                            if debug_backoff; fprintf('[%d us @ %d] channel is idle, current state: STATE_CSMA.\n',(i-1)*10, j); end
                            nodeState(j)=STATE_CSMA;
                            delay(j)=delay(j)-1;
                            if delay(j)==0
        %                         temp = sprintf('id %d csma1', j);
        %                         disp(temp);
                                cnt=cnt+1;
                                if cnt>1
                                    if debug_collision; fprintf('[%d us @ %d] collision!!!!!!!!!!!!!!\n',(i-1)*10, j); end
                                    nbCollision(j)=nbCollision(j)+1;
                                    %nbTransmission(j)=nbTransmission(j)+1;
                                    if stateintime == TIMESTATE_DATA
                                        tempbusyFor(j)=LDATA+LACKTIMEOUT;
                                    elseif stateintime == TIMESTATE_ATIM
                                        tempbusyFor(j)=LATIM+LACKTIMEOUT;
                                    end
                                    trfail(j)=1;
                                elseif cnt==1
                                    rd=j;
                                end
    %                             fprintf('[%d] %d will begin to transmit.',i,j);
                                nodeState(j)=STATE_TX;
                            end
                            %disp('channel idle');
                        else % channle not idle
                            if debug_backoff; fprintf('[%d us @ %d] channel is busy, current state: STATE_CSMA.\n',(i-1)*10, j); end
                            oneshot(j)=0;
                            DIFScnt(j)=0;
                            nodeState(j)=STATE_DIFS;
                        end
                    case{STATE_TX}
                        tempbusyFor(j) = busyFor(j);% bug!!! fixed on 2/13/15 (almost 3 years!)
                        tempbusyFor(j) = tempbusyFor(j) - 1;% bug!!! fixed on 2/13/15 (almost 3 years!)
                        if tempbusyFor(j)==0
    %                         tempbusyFor(j)=0;
                            nbTransmission(j)=nbTransmission(j)+1;
                            if trfail(j)==0% so this means this transmission is collision free
                                successTX(j)=1; %and also it is finished, so it is successful.
                                if debug_tx; fprintf('[%d us @ %d] send successfully!\n',(i-1)*10, j); end
                                
                                if i - StartTimeofCurrentCycle > BEACON_WINDOW + ATIM_WINDOW && i - StartTimeofCurrentCycle < TOTAL_WINDOW
                                    if debug_psm; fprintf('[%d us @ %d] DATA sends OK, send new ATIM+DATA in next round... awaketime: %f += %f\n', (i-1)*10, j, Awaketime(j), i - StartTimeofCurrentCycle); end
%                                     nbDataTransmission(j) = nbDataTransmission(j)+1;
                                    Awaketime(j) = Awaketime(j) + (i - StartTimeofCurrentCycle);
                                end

                            else
                                if debug_tx; fprintf('[%d us @ %d] send unsuccessfully.\n',(i-1)*10, j); end
                            end
    %                         if SATURATION_ENABLED == 1
                                if successTX(j)==1
                                    nodeState(j)=STATE_WAIT;
                                else
                                    nodeState(j)=STATE_INIT;
                                end
    %                         else
    %                             nodeState(j)=STATE_WAIT;
    %                         end
    %                         fprintf('[%d us] node %d transmit done.\n',(i-1)*10,j);
                        else
                            if debug_tx; fprintf('[%d us @ %d] is transmitting! busyFor: %d, tempbusyFor: %d\n', ...
                                    (i-1)*10, j, busyFor(j), tempbusyFor(j)); end
                        end
                end

                if DIFScnt(j)==DIFS
                    DIFScnt(j)=0;
    %                 if oneshot(j)==1
    %     %                 fprintf('id %d oneshot', j);
    % 
    %                     cnt=cnt+1;
    %                     if cnt>1
    %                         fprintf('[%d us] collision!!!!!!!!!!!!!!: node %d\n',(i-1)*10, j);
    %                         nbCollision(j)=nbCollision(j)+1;
    %                         %nbTransmission(j)=nbTransmission(j)+1;
    %                         tempbusyFor(j)=LDATA+LACKTIMEOUT;
    %                         trfail(j)=1;
    %                     elseif cnt==1
    %                         rd=j;
    %                     end
    %                     nodeState(j)=STATE_TX;
    % %                     fprintf('[%d us] node %d will begin to transmit(!).\n',(i-1)*10,j);
    % 
    %                 else
                        if delay(j)==0
                            if trfail(j)==1
                                CW(j)=min(aMaxBE,CW(j)+1);% aMinBE < BE < aMaxBE
                                delay(j)=int32(rand*2^(CW(j))); % must add 1, 
                                if debug_backoff; fprintf('[%d us @ %d] backoff (%d) for %d\n', (i-1)*10, j, CW(j), delay(j)); end
                            else
                                CW(j)=aMinBE; % aMinBE < BE < aMaxBE
                                delay(j)=int32(rand*2^(CW(j))); % must add 1,
                                if debug_backoff; fprintf('[%d us @ %d] backoff (%d) for %d\n', (i-1)*10, j, CW(j), delay(j)); end
                            end
    %                         fprintf('[%d] %d begin to delay for %d.',i,j,delay(j));

        %                     fprintf('id %d delay %d', j, delay(j));

                            if delay(j)==0
                                cnt=cnt+1;
        %                         fprintf('id %d csma0', j);

                                if cnt>1
                                    if debug_collision; fprintf('[%d us @ %d] collision!!!!!!!!!!!!!!\n',(i-1)*10, j); end
                                    nbCollision(j)=nbCollision(j)+1;
                                    %nbTransmission(j)=nbTransmission(j)+1;
                                    if stateintime == TIMESTATE_DATA
                                        tempbusyFor(j)=LDATA+LACKTIMEOUT;
                                    elseif stateintime == TIMESTATE_ATIM
                                        tempbusyFor(j)=LATIM+LACKTIMEOUT;
                                    end
                                    trfail(j)=1;
                                elseif cnt==1
                                    rd=j;
                                end
                                nodeState(j)=STATE_TX;
    %                             temp=sprintf('[%d] %d will begin to transmit.',i,j);
    %                             disp(temp);
                            else
                                if sum(busyFor(1:N))==0 % channel idle
                                    nodeState(j)=STATE_CSMA;
                                    %delay(j)=delay(j)-1;
                                    %disp('channel idle');
                                else % channle not idle
                                    oneshot(j)=0;
                                    DIFScnt(j)=0;
                                    nodeState(j)=STATE_DIFS;
                                end
                            end

                        else
    %                         temp=sprintf('[%d] %d continue to delay for %d.',i,j,delay(j));
    %                         disp(temp);
                            if sum(busyFor(1:N))==0 % channel idle
                                nodeState(j)=STATE_CSMA;
                                %delay(j)=delay(j)-1;
                                %disp('channel idle');
                            else % channle not idle
                                oneshot(j)=0;
                                DIFScnt(j)=0;
                                nodeState(j)=STATE_DIFS;
                            end      
                        end
    %                 end
                end
            end

            if cnt>1
                nbCollision(rd)=nbCollision(rd)+1;
                if debug_collision; fprintf('[%d us @ %d] collision!!!!!!!!!!!!!!\n', (i-1)*10, rd); end
                %nbTransmission(rd)=nbTransmission(rd)+1;
                busyFor=tempbusyFor;
                if stateintime == TIMESTATE_DATA
                    busyFor(rd)=LDATA+LACKTIMEOUT;
                elseif stateintime == TIMESTATE_ATIM
                    busyFor(rd)=LATIM+LACKTIMEOUT;
                end
    %             busyFor(rd)=LDATA+LACKTIMEOUT;
                tempbusyFor=busyFor;% bug!!! fixed on 2/13/15 (almost 3 years!)
                trfail(rd)=1;
            elseif cnt==1
                %nbTransmission(rd)=nbTransmission(rd)+1;
                if stateintime == TIMESTATE_DATA
                    busyFor(rd)=LDATA+LACK;
                elseif stateintime == TIMESTATE_ATIM
                    busyFor(rd)=LATIM+LACK;
                end
    %             busyFor(rd)=LDATA+LACK;
                tempbusyFor=busyFor;% bug!!! fixed on 2/13/15 (almost 3 years!)
                trfail(rd)=0;
            end
            cnt=0;

            for j=1:N
                busyFor(j)=tempbusyFor(j);% bug!!! fixed on 2/13/15 (almost 3 years!)
            end

            if DYNAMIC_PLOT == 1
                if mod(i,UPDATE_INTERVAL)==0
    %                 disp('update drawing!!!!!');
                    m=[m n];
    %                 disp(size(m));
                    if size(m,2)==SHOW_INTERVAL/UPDATE_INTERVAL+2
                        m=m(:,2:SHOW_INTERVAL/UPDATE_INTERVAL+1+1);
    %                     disp(size(m,1));
    %                     disp(size(m,2));
    %                     return;
                    end

                    for j=1:N
                        set(p(j),'XData',t,'YData',m(j,:));
                    end
                    drawnow;

                end
            end


            global_sync = global_sync + 1;

        end
        diary on;
        %disp(sum(nbTransmission)/N);
    %     fprintf('N: %d\n', N);
%         disp(Awaketime);
%         disp(nbDataTransmission);
        fprintf('average awake time: %f ms\n', sum(Awaketime)/N/(T/TOTAL_WINDOW)/100);
        fprintf('# pernode data TXs: %f\n', sum(nbDataTransmission)/N);
%         disp(TXtime./nbDataTransmission);
        fprintf('average perpacket TX time: %f ms\n', sum((TXtime./nbDataTransmission))/N/100);
        fprintf('queue empty prob.: %f\n', sum(queue_empty_count)/N/T);
        rho = traffic_rate/(1000/(sum((TXtime./nbDataTransmission))/N/100));
        fprintf('average queue size: %f, theoretical: %f\n', sum(queue_size_count)/N/T, rho+(rho^2)/(2*(1-rho)));
        averageTrans(N)=sum(nbDataTransmission)/N*PAYLOAD*8/(T*10/1000000)/1000000;
        fprintf('average throughput: %f Mbps\n', averageTrans(N));
        totalTrans(N)=sum(nbDataTransmission)*PAYLOAD*8/(T*10/1000000)/1000000;
        fprintf('total throughput: %f Mbps\n', totalTrans(N));

    %     averageTrans(N)=sum(nbTransmission-nbCollision)/N*PAYLOAD*8/(T*10/1000000)/1000000;
    %     disp(['average throughput: ', num2str(averageTrans(N))]);
    %     totalTrans(N)=sum(nbTransmission-nbCollision)*PAYLOAD*8/(T*10/1000000)/1000000;
    %     disp(['total throughput: ', num2str(totalTrans(N))]);
    %     fprintf('sum(nbTransmission): %g, sum(nbTransmission)/%d: %g\n', sum(nbTransmission), N, sum(nbTransmission)/N);
    %     fprintf('sum(nbCollision): %g sum(nbCollision)/%d: %g\n', sum(nbCollision), N, sum(nbCollision)/N);
    %     
    %     for j=1:N
    %         collisions(j) = nbCollision(j)/nbTransmission(j);
    %     end
    %     fprintf('collisions rate: %g\n', sum(collisions)/N);

    end
        fprintf('---- Running time=%g s ----\n', etime(clock, tstart));
    diary off;

    % figure;
    % % plot(averageTrans,'DisplayName','averageTrans','YDataSource','averageTrans');hold all;
    % % plot(totalTrans,'DisplayName','totalTrans','YDataSource','averageTrans');grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,50,0,50]);hold all;figure(gcf);
    % plot(collisions,'DisplayName','averageTrans','YDataSource','averageTrans');hold all;
end