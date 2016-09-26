%function[]=wifi()
%time slot=10us


%State Vectors
%N=1;
aMinBE=5;
aMaxBE=10;
LDATA=7;
LMAC=2;
LMACTIMEOUT=4;
T=100000;
%int32(rand.*(2.^(BE(1:N))));
%for i=1:N
%    delay(i)=int32(rand*2^(BE(i)));
%end

%nbFailure(1:N)=0;
%nbCCA(1:N)=0;

DIFS=3;
STATE_INIT=0;
STATE_DIFS=1;
STATE_CSMA=2;
STATE_TRDELAY=3;
STATE_SOFTWARE=4; % for software delay
%STATE_RETR=4;
%disp(delay(1:N));
cnt=0;
rd=0;

Trials=6;
averageCCA(1:Trials)=0;
averageBusy1(1:Trials)=0;
averageBusy2(1:Trials)=0;
averageTrans(1:Trials)=0;
averageSucTrans(1:Trials)=0;

for N=1:Trials
    %disp(N);
    
    delay(1:Trials)=0;
    CW(1:Trials)=aMinBE;% aMinBE < BE < aMaxBE
    busyFor(1:Trials)=0; %for transmission delay
    tempbusyFor(1:Trials)=0; %for transmission delay
    nbTransmission(1:Trials)=0;
    nbCollision(1:Trials)=0;
    DIFScnt(1:Trials)=0; %for DIFS delay
    oneshot(1:Trials)=1; %for oneshot of transmission, 1 means the first try, 0 means the first try fails.
    nodeState(1:Trials)=0; %for FSM.
    trfail(1:Trials)=0; %for tr fails.
    softdelay(1:Trials)=10;
    
    for i=1:T
        for j=1:N
            if nodeState(j)==STATE_INIT
               DIFScnt(j)=0;
                if trfail(j)==0
                    oneshot(j)=1;
                else
                    oneshot(j)=0;
                end
                nodeState(j)=STATE_DIFS;
            elseif nodeState(j)==STATE_SOFTWARE
                softdelay(j)=softdelay(j)-1;
                if softdelay(j)==0
                    nodeState(j)=STATE_INIT;
                    softdelay(j)=10;
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
                                nbCollision(j)=nbCollision(j)+1;
                                %nbTransmission(j)=nbTransmission(j)+1;
                                tempbusyFor(j)=LDATA+LMACTIMEOUT;
                                trfail(j)=1;
                            elseif cnt==1
                                rd=j;
                            end
                            temp=sprintf('[%d] %d will begin to transmit.',i,j);
                            disp(temp);
                            nodeState(j)=STATE_TRDELAY;
                        end
                        %disp('channel idle');
                    else % channle not idle
                        oneshot(j)=0;
                        DIFScnt(j)=0;
                        nodeState(j)=STATE_DIFS;
                    end
                case{STATE_TRDELAY}
                    busyFor(j)=busyFor(j)-1;
                    if busyFor(j)==0
                        tempbusyFor(j)=0;
                        nbTransmission(j)=nbTransmission(j)+1;
                        nodeState(j)=STATE_SOFTWARE;
                        temp=sprintf('[%d] %d transmit done.',i,j);
                        disp(temp);
                    end
            end

            if DIFScnt(j)==DIFS
                DIFScnt(j)=0;
                if oneshot(j)==1
    %                 temp = sprintf('id %d oneshot', j);
    %                 disp(temp);
    %                     temp = sprintf('fuck %d', j);
    %                     temp = sprintf('you %s %d', temp, j);
    %                     disp(temp);
                    cnt=cnt+1;
                    if cnt>1
                        nbCollision(j)=nbCollision(j)+1;
                        %nbTransmission(j)=nbTransmission(j)+1;
                        tempbusyFor(j)=LDATA+LMACTIMEOUT;
                        trfail(j)=1;
                    elseif cnt==1
                        rd=j;
                    end
                    nodeState(j)=STATE_TRDELAY;
                    temp=sprintf('[%d] %d will begin to transmit(!).',i,j);
                    disp(temp);
                else
                    if delay(j)==0
                        if trfail(j)==1
                            CW(j)=min(aMaxBE,CW(j)+1);% aMinBE < BE < aMaxBE
                            delay(j)=int32(rand*2^(CW(j)))+1; % must add 1, 
                        else
                            CW(j)=aMinBE; % aMinBE < BE < aMaxBE
                            delay(j)=int32(rand*2^(CW(j)))+1; % must add 1,
                        end
                        temp=sprintf('[%d] %d begin to delay for %d.',i,j,delay(j));
                        disp(temp);

    %                     temp = sprintf('id %d delay %d', j, delay(j));
    %                     disp(temp);

                        if delay(j)==0
                            cnt=cnt+1;
    %                         temp = sprintf('id %d csma0', j);
    %                         disp(temp);
                            if cnt>1
                                nbCollision(j)=nbCollision(j)+1;
                                %nbTransmission(j)=nbTransmission(j)+1;
                                tempbusyFor(j)=LDATA+LMACTIMEOUT;
                                trfail(j)=1;
                            elseif cnt==1
                                rd=j;
                            end
                            nodeState(j)=STATE_TRDELAY;
                            temp=sprintf('[%d] %d will begin to transmit.',i,j);
                            disp(temp);
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
                        temp=sprintf('[%d] %d continue to delay for %d.',i,j,delay(j));
                        disp(temp);
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
            %nbTransmission(rd)=nbTransmission(rd)+1;
            busyFor=tempbusyFor;
            busyFor(rd)=LDATA+LMACTIMEOUT;
            trfail(rd)=1;
        elseif cnt==1
            %nbTransmission(rd)=nbTransmission(rd)+1;
            busyFor(rd)=LDATA+LMAC;
            trfail(rd)=0;
        end
        cnt=0;
    end
    
    %disp(sum(nbTransmission)/N);
    averageTrans(N)=sum(nbTransmission-nbCollision)/N*512*8/(T*10/1000000)/1000000;
    disp(averageTrans(N));
    
end


plot(averageTrans,'DisplayName','averageTrans','YDataSource','averageTrans');grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,50,0,50]);figure(gcf);
