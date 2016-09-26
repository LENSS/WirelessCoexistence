%N=1;

Trials=20;
% averageCCA(1:Trials)=0;
% averageBusy1(1:Trials)=0;
% averageBusy2(1:Trials)=0;
% averageTrans(1:Trials)=0;
% averageSucTrans(1:Trials)=0;
averageWifiTrans(1:Trials)=0;
averageBmacTrans(1:Trials)=0;

for N=1:Trials

%time slot=10us
    wifiTRcount(1:Trials)=0;
    bmacTRcount(1:Trials)=0;
    wifiCollision(1:Trials)=0;
    bmacCollision(1:Trials)=0;

    aMinBE=5;
    aMaxBE=10;
    CW(1:Trials)=aMinBE;% aMinBE < BE < aMaxBE
    LDATA=7;
    LMAC=2;
    LMACTIMEOUT=4;
    delay(1:Trials)=0;
    %int32(rand.*(2.^(BE(1:N))));
    %for i=1:N
    %    delay(i)=int32(rand*2^(BE(i)));
    %end

    busyFor(1:Trials)=0; %for transmission delay
    tempbusyFor(1:Trials)=0; %for transmission delay
    T=1000000;
    nbTransmission(1:Trials)=0;
    nbCollision(1:Trials)=0;
    %nbFailure(1:N)=0;

    DIFScnt(1:Trials)=0; %for DIFS delay
    oneshot(1:Trials)=1; %for oneshot of transmission, 1 means the first try, 0 means the first try fails.
    nodeState(1:Trials)=0; %for FSM.
    trfail(1:Trials)=0; %for tr fails.
    DIFS=3;
    STATE_CHOOSE=0; % WIFI & BMAC
    STATE_INIT=1;
    STATE_DIFS=2;
    STATE_CSMA=3;
    STATE_TRDELAY=4;
    STATE_SOFTWARE=5; % for software delay
    %STATE_RETR=4;
    %disp(delay(1:N));
    cnt=0;
    rd=0;
    c=0;
    nodeType(1:Trials)=0; %for WIFI/BMAC.
    NODE_WIFI=0;
    NODE_BMAC=1;
    softdelay(1:Trials)=10;



    % aMinBEb=3;
    % aMaxBEb=5;
    % macMaxCSMABackoffs=5;
    %NBb(1:N)=0;% <macMaxCSMABackoffs
    CWb(1:Trials)=2;% 
    %BEb(1:N)=aMinBE;% aMinBE < BE < aMaxBE
    nbCCA(1:Trials)=0;
    LDATAb=128*3;
    softwareDelay(1:Trials)=0;
    wificount=0;
    bmaccount=0;
    if N==1
        wificount=1;
    elseif N==2
        wificount=2;
    else
        wificount=2;
        bmaccount=N-2;
    end


    for i=1:T
        for j=1:N
            if nodeState(j)==STATE_CHOOSE
                p=binornd(1,0.8); % binornd(1,x) x is the probability of bmac behavier.
                %if p==0 % wifi 0.2
                if j==1 || j==2
                    nodeType(j)=NODE_WIFI;
                    nodeState(j)=STATE_INIT;
                else %bmac 0.8
                    nodeType(j)=NODE_BMAC;
                    nodeState(j)=STATE_CSMA;
                    delay(j)=int32(10+rand*310)*3;
                end
            end
            %disp('1');
            if nodeType(j)==NODE_WIFI
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
                                cnt=cnt+1;
                                if cnt>1
                                    nbCollision(j)=nbCollision(j)+1;
                                    nbTransmission(j)=nbTransmission(j)+1;
                                    wifiCollision(j)=wifiCollision(j)+1;
                                    wifiTRcount(j)=wifiTRcount(j)+1;
                                    tempbusyFor(j)=LDATA+LMACTIMEOUT;
                                    trfail(j)=1;
                                elseif cnt==1
                                    rd=j;
                                end
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
                            %delay(j)=0;
                            nodeState(j)=STATE_SOFTWARE;
                        end
                end

                if DIFScnt(j)==DIFS
                    DIFScnt(j)=0;
                    if oneshot(j)==1
        %                     temp = sprintf('fuck %d', j);
        %                     temp = sprintf('you %s %d', temp, j);
        %                     disp(temp);
                        cnt=cnt+1;
                        if cnt>1
                            nbCollision(j)=nbCollision(j)+1;
                            nbTransmission(j)=nbTransmission(j)+1;
                            wifiCollision(j)=wifiCollision(j)+1;
                            wifiTRcount(j)=wifiTRcount(j)+1;
                            tempbusyFor(j)=LDATA+LMACTIMEOUT;
                            trfail(j)=1;
                        elseif cnt==1
                            rd=j;
                        end
                        nodeState(j)=STATE_TRDELAY;
                    else
                        if delay(j)==0
                            if trfail(j)==1
                                CW(j)=min(aMaxBE,CW(j)+1);% aMinBE < BE < aMaxBE
                                delay(j)=int32(rand*2^(CW(j)))+1; % must add 1, 
                            else
                                CW(j)=aMinBE; % aMinBE < BE < aMaxBE
                                delay(j)=int32(rand*2^(CW(j)))+1; % must add 1,
                            end

                            if delay(j)==0
                                cnt=cnt+1;
                                if cnt>1
                                    nbCollision(j)=nbCollision(j)+1;
                                    nbTransmission(j)=nbTransmission(j)+1;
                                    wifiCollision(j)=wifiCollision(j)+1;
                                    wifiTRcount(j)=wifiTRcount(j)+1;
                                    tempbusyFor(j)=LDATA+LMACTIMEOUT;
                                    trfail(j)=1;
                                elseif cnt==1
                                    rd=j;
                                end
                                nodeState(j)=STATE_TRDELAY;
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
        end

        for j=1:N
            if nodeType(j)==NODE_BMAC
                if softwareDelay(j)~=0
                    continue;
                end
                if delay(j)==0% && lock(j)==0
                    nbCCA(j)=nbCCA(j)+1;
                end
            end
        end

        if sum(busyFor(1:N))==0 % channel idle
            %disp('channel idle');
            for j=1:N
                if nodeType(j)==NODE_BMAC
                    if softwareDelay(j)~=0
                        continue;
                    end
                    if delay(j)==0 || delay(j)==-3
                        %disp('cca');
                        CWb(j)=CWb(j)-1;
                        if CWb(j)==0
                            cnt=cnt+1;
                            if cnt>1
                                nbCollision(j)=nbCollision(j)+1;
                                nbTransmission(j)=nbTransmission(j)+1;
                                bmacCollision(j)=bmacCollision(j)+1;
                                bmacTRcount(j)=bmacTRcount(j)+1;
                                tempbusyFor(j)=LDATAb;
                            elseif cnt==1
                                rd=j;
                            end
                        end
                    end
                end
            end
%             if cnt>1
%                 nbCollision(rd)=nbCollision(rd)+1;
%                 nbTransmission(rd)=nbTransmission(rd)+1;
%                 busyFor(rd)=busyFor(rd)+length;
%             elseif cnt==1
%                 nbTransmission(rd)=nbTransmission(rd)+1;
%                 busyFor(rd)=busyFor(rd)+length;
%             end
        else % channle not idle
            %disp('channel not idle');
            for j=1:N
                if nodeType(j)==NODE_BMAC
                    if softwareDelay(j)~=0
                        continue;
                    end
                   if busyFor(j)>0
                        busyFor(j)=busyFor(j)-1;
                        if busyFor(j)==0
                            %disp('sendone');
                            %NBb(j)=0;% <macMaxCSMABackoffs
                            CWb(j)=2;% 
                            %BEb(j)=aMinBEb;% aMinBE < BE < aMaxBE
                            delay(j)=int32(10+rand*310)*3+1; % must add 1, 
                            nodeState(j)=STATE_CHOOSE;
                            %because this time slot has already used to do busyFor
                            softwareDelay(j)=210*3; % 6.5ms
                        end
                    end
                end
           end

            for j=1:N
                if nodeType(j)==NODE_BMAC
                    if softwareDelay(j)~=0
                        continue;
                    end
                    if delay(j)==0 || delay(j)==-3
                        %disp('backoff');
        %                 NB(j)=NB(j)+1;% <macMaxCSMABackoffs
                        CWb(j)=2;% 
        %                 BE(j)=min(aMaxBE,BE(j)+1);% aMinBE < BE < aMaxBE
                        %BEb(j)=aMaxBE;% aMinBE < BE < aMaxBE
                        delay(j)=int32(10+rand*70)*3+1; % must add 1, 
                        %because this time slot has already used to do update
        %                if NB(j)==macMaxCSMABackoffs
        %                    nbFailure(j)=nbFailure(j)+1;
        %                    NB(j)=0;% <macMaxCSMABackoffs
        %                    CW(j)=2;% 
        %                    BE(j)=aMinBE;% aMinBE < BE < aMaxBE
        %                    delay(j)=int32(rand*2^(BE(j)))+1; % must add 1, 
        %                    %because this time slot has already used to do reset
        %                end
                    end
                end
            end

        end
        %cnt=0;

        for j=1:N
            if nodeType(j)==NODE_BMAC
                if softwareDelay(j)~=0
                    softwareDelay(j)=softwareDelay(j)-1;
                    continue;
                end
                if delay(j)~=-4
                    delay(j)=delay(j)-1;
                end
            end
        end

        if cnt>1
            nbCollision(rd)=nbCollision(rd)+1;
            nbTransmission(rd)=nbTransmission(rd)+1;
            busyFor=tempbusyFor;
            if nodeType(rd)==NODE_WIFI
                wifiCollision(rd)=wifiCollision(rd)+1;
                wifiTRcount(rd)=wifiTRcount(rd)+1;
                busyFor(rd)=LDATA+LMACTIMEOUT;
            else
                bmacCollision(rd)=bmacCollision(rd)+1;
                bmacTRcount(rd)=bmacTRcount(rd)+1;
                busyFor(rd)=LDATAb;
            end
            trfail(rd)=1;
        elseif cnt==1
            nbTransmission(rd)=nbTransmission(rd)+1;
            if nodeType(rd)==NODE_WIFI
                wifiTRcount(rd)=wifiTRcount(rd)+1;
                busyFor(rd)=LDATA+LMAC;
            else
                bmacTRcount(rd)=bmacTRcount(rd)+1;
                busyFor(rd)=LDATAb;
            end
           trfail(rd)=0;
        end

        cnt=0;    

    end
    
    disp(wificount);disp(bmaccount);
    averageWifiTrans(N)=sum(wifiTRcount-wifiCollision)/wificount*512*8/(T*10/1000000)/1000000;
    disp(averageWifiTrans(N));
    averageBmacTrans(N)=sum(bmacTRcount-bmacCollision)/bmaccount*8*100/(T*10/1000000)/1000;
    disp(averageBmacTrans(N));
    
end

