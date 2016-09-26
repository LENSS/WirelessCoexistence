
clear all;

rng('default');
rng(0);

% Create a normally distributed (mu: 5, sigma: 3) random data set
% x = normrnd(5, 3, 1e4, 1);
x = exprnd(2, 1e5, 1);
y = zeros(1e4, 1);
for i=1:1e4
    y(i) = ceil(x(i*10)*1e5)/1e5;
end

% Compute and plot results. The results are sorted by "Bayesian information
% criterion".
[D, PD] = allfitdist(y, 'PDF');
D(1)

return;


%State Vectors
PAYLOAD = 1000;%bytes
HEDEAR = 40;%bytes, mac header, tcp/udp header and ip header
PREAMBLE = 4;%40 us
aMinBE=5;
aMaxBE=10;

LDATA=4+ceil(HEDEAR*8/54/10)+ceil(PAYLOAD*8/54/10);%ceil(846*8/54/10)
LACK=2;
LACKTIMEOUT=8;


T=100000;

%nbFailure(1:N)=0;
%nbCCA(1:N)=0;

DIFS=3;
STATE_IDLE=0; % for waiting for a packet
STATE_SOFT=1; % for software delay
STATE_INIT=2;
STATE_DIFS=3;
STATE_CSMA=4;
STATE_TX=5;

cnt=0;
rd=0;

queue_size = 50;
SOFT_DELAY = 1;
SATURATION_ENABLED = 0; %0: unsaturated
if SATURATION_ENABLED == 1
    DYNAMIC_PLOT = 0;
else
    DYNAMIC_PLOT = 1;
end
UPDATE_INTERVAL = 1000; %for ploting, not for traffic!

traffic_rate = 200; %packets/s

TR_POSSION = 1;
TR_PERIOD = 2;
TRAFFIC_TYPE = TR_POSSION;%TR_PERIOD;%


Trials=20;%number of nodes


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
    nbCollision(1:Trials)=0;
    DIFScnt(1:Trials)=0; %for DIFS delay
    oneshot(1:Trials)=0; %for oneshot of transmission, 1 means the first try, 0 means the first try fails.
    if SATURATION_ENABLED == 1
        nodeState(1:Trials)=STATE_SOFT; %for FSM 
    else
        nodeState(1:Trials)=STATE_IDLE; %for FSM 
    end
    trfail(1:Trials)=0; %for tr fails.
    softdelay(1:Trials)=SOFT_DELAY;
    interarrival(1:Trials)=-1;
    
    
    
    mac_queue(1:Trials)=queue_size;
    queue_usage(1:Trials)=0;
    
    if DYNAMIC_PLOT == 1

        t=[0];
        m=[];
        for i=1:N
            m=[m;0];
        end

        p = plot(t,m,...
           'EraseMode','background','MarkerSize',5);hold on;
        axis([1 T/UPDATE_INTERVAL -5 queue_size+5]);
        grid on;

    end    
    
    
    tstart = clock;

    for i=1:T
        
        if DYNAMIC_PLOT == 1
                
            %draw figure to show the queue usage
            if mod(i,UPDATE_INTERVAL)==0
                t=[t i/UPDATE_INTERVAL];                   %Matrix 1*(i+1)
                n=[];
            end
        end
        
        for j=1:N
            
            if DYNAMIC_PLOT == 1
                if mod(i,UPDATE_INTERVAL)==0
                    n=[n;queue_size-mac_queue(j)];
                end
            end
            
            %Traffic generator
            %===============================================================================
            if SATURATION_ENABLED == 0

                if interarrival(j)==-1
                    while(1)
                        
                        if TRAFFIC_TYPE == TR_POSSION
                            interarrival(j) = ceil(exprnd(1/traffic_rate)*1e5);
                        elseif TRAFFIC_TYPE == TR_PERIOD
                            interarrival(j) = ceil((1/traffic_rate)*1e5);
                        end
                        
                        if interarrival(j)==0
                            %put data into mac_queue
                            if mac_queue(j)==0
                                fprintf('[%d us] !!!!!!no space in the mac_queue(%d), drop it!!!!!\n', (i-1)*10, j);
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
                        fprintf('[%d us] no space in the mac_queue(%d), drop it!!!!!\n', (i-1)*10, j);
                    else
                        mac_queue(j) = mac_queue(j) - 1;
%                         fprintf('[%d us] put data into the mac_queue(%d)\n', (i-1)*10, j);
                    end

                    while(1)
                        interarrival(j) = ceil(exprnd(1/traffic_rate)*1e5);
                        if interarrival(j)==0
                        %put data into mac_queue
                            if mac_queue(j)==0
                                fprintf('[%d us] no space in the mac_queue(%d), drop it!!!!!\n', (i-1)*10, j);
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
            end
            %===============================================================================

            
%             continue;
            
            
            if nodeState(j)==STATE_IDLE
                %check if mac_queue has data

                if mac_queue(j)~=queue_size
                    %fetch data from mac_queue
                    mac_queue(j) = mac_queue(j) + 1;
%                     fprintf('[%d us] get a data in mac_queue(%d)\n', (i-1)*10, j);
                    nodeState(j)=STATE_SOFT;
                else
                end
                
            elseif nodeState(j)==STATE_INIT
               DIFScnt(j)=0;
                if trfail(j)==0
                    oneshot(j)=0;
                else
                    oneshot(j)=0;
                end
                nodeState(j)=STATE_DIFS;
                continue;
            elseif nodeState(j)==STATE_SOFT
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
                        DIFScnt(j)=DIFScnt(j)+1;
                    else % channle not idle
                        oneshot(j)=0;
                        DIFScnt(j)=0;
                    end
                case{STATE_CSMA}
                    if sum(busyFor(1:N))==0 % channel idle
                        nodeState(j)=STATE_CSMA;
                        delay(j)=delay(j)-1;
                        if delay(j)==0
    %                         temp = sprintf('id %d csma1', j);
    %                         disp(temp);
                            cnt=cnt+1;
                            if cnt>1
                                fprintf('[%d us] 1 collision!!!!!!!!!!!!!!: node %d\n',(i-1)*10, j);
                                nbCollision(j)=nbCollision(j)+1;
                                %nbTransmission(j)=nbTransmission(j)+1;
                                tempbusyFor(j)=LDATA+LACKTIMEOUT;
                                trfail(j)=1;
                            elseif cnt==1
                                rd=j;
                            end
%                             fprintf('[%d] %d will begin to transmit.',i,j);
                            nodeState(j)=STATE_TX;
                        end
                        %disp('channel idle');
                    else % channle not idle
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
                        if SATURATION_ENABLED == 1
                            nodeState(j)=STATE_SOFT;
                        else
                            nodeState(j)=STATE_IDLE;
                        end
%                         fprintf('[%d us] node %d transmit done.\n',(i-1)*10,j);

                    end
            end

            if DIFScnt(j)==DIFS
                DIFScnt(j)=0;
                if oneshot(j)==1
    %                 fprintf('id %d oneshot', j);

                    cnt=cnt+1;
                    if cnt>1
                        fprintf('[%d us] 2 collision!!!!!!!!!!!!!!: node %d\n',(i-1)*10, j);
                        nbCollision(j)=nbCollision(j)+1;
                        %nbTransmission(j)=nbTransmission(j)+1;
                        tempbusyFor(j)=LDATA+LACKTIMEOUT;
                        trfail(j)=1;
                    elseif cnt==1
                        rd=j;
                    end
                    nodeState(j)=STATE_TX;
%                     fprintf('[%d us] node %d will begin to transmit(!).\n',(i-1)*10,j);

                else
                    if delay(j)==0
                        if trfail(j)==1
                            CW(j)=min(aMaxBE,CW(j)+1);% aMinBE < BE < aMaxBE
                            delay(j)=int32(rand*2^(CW(j))); % must add 1, 
%                             fprintf('[%d us] node %d backoff (%d) for %d\n', (i-1)*10, j, CW(j), delay(j)); 
                        else
                            CW(j)=aMinBE; % aMinBE < BE < aMaxBE
                            delay(j)=int32(rand*2^(CW(j))); % must add 1,
%                             fprintf('[%d us] node %d backoff (%d) for %d\n', (i-1)*10, j, CW(j), delay(j)); 
                        end
%                         fprintf('[%d] %d begin to delay for %d.',i,j,delay(j));

    %                     fprintf('id %d delay %d', j, delay(j));

                        if delay(j)==0
                            cnt=cnt+1;
    %                         fprintf('id %d csma0', j);

                            if cnt>1
                                fprintf('[%d us] 3 collision!!!!!!!!!!!!!!: node %d\n',(i-1)*10, j);
                                nbCollision(j)=nbCollision(j)+1;
                                %nbTransmission(j)=nbTransmission(j)+1;
                                tempbusyFor(j)=LDATA+LACKTIMEOUT;
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
                end
            end
        end
        if cnt>1
            nbCollision(rd)=nbCollision(rd)+1;
            fprintf('[%d us] 4 collision!!!!!!!!!!!!!!: node %d\n', (i-1)*10,rd);
            %nbTransmission(rd)=nbTransmission(rd)+1;
            busyFor=tempbusyFor;
            busyFor(rd)=LDATA+LACKTIMEOUT;
            tempbusyFor=busyFor;% bug!!! fixed on 2/13/15 (almost 3 years!)
            trfail(rd)=1;
        elseif cnt==1
            %nbTransmission(rd)=nbTransmission(rd)+1;
%             fprintf('[%d us] successful tx!!!!!!!!!!!!!!: node %d\n', (i-1)*10,rd);
            busyFor(rd)=LDATA+LACK;
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

                for j=1:N
                    set(p(j),'XData',t,'YData',m(j,:))
                end

                drawnow;
            end
        end
        
        
        
    end
    
    %disp(sum(nbTransmission)/N);
    fprintf('N: %d\n', N);
    averageTrans(N)=sum(nbTransmission-nbCollision)/N*PAYLOAD*8/(T*10/1000000)/1000000;
    disp(['average throughput: ', num2str(averageTrans(N))]);
    totalTrans(N)=sum(nbTransmission-nbCollision)*PAYLOAD*8/(T*10/1000000)/1000000;
    disp(['total throughput: ', num2str(totalTrans(N))]);
    fprintf('sum(nbTransmission)/N: %g\n', sum(nbTransmission)/N);
    fprintf('sum(nbCollision)/N: %g\n', sum(nbCollision)/N);
    collisions(N) = sum(nbCollision)/sum(nbTransmission);
    disp(['collisions rate: ', num2str(collisions(N))]);    
end
    fprintf('---- Running time=%g s ----\n', etime(clock, tstart));

figure;
% plot(averageTrans,'DisplayName','averageTrans','YDataSource','averageTrans');hold all;
% plot(totalTrans,'DisplayName','totalTrans','YDataSource','averageTrans');grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,50,0,50]);hold all;figure(gcf);
plot(collisions,'DisplayName','averageTrans','YDataSource','averageTrans');hold all;
