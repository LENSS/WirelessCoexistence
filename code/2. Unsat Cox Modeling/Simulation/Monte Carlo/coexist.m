
function ret = coexist(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b)
%N=1;
% clear all;

% [T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b]=deal(1000000, 1, 1, 40, 40, 16, 1024, 1500, 310, 70, 108);


STEP_SIZE=5;
TotalTrials = 200;
Statistics = 0;
if ~(BMAC_END==0 && BMAC_START==0)
    Statistics = Statistics + (BMAC_END-BMAC_START)/STEP_SIZE+1;
end

if ~(WIFI_END==0 && WIFI_START==0)
    Statistics = Statistics + (WIFI_END-WIFI_START)/STEP_SIZE+1;
end


% averageCCA(1:Trials)=0;
% averageBusy1(1:Trials)=0;
% averageBusy2(1:Trials)=0;
% averageTrans(1:Trials)=0;
% averageSucTrans(1:Trials)=0;
% averageWifiTrans(1:Statistics)=0;
% averageBmacTrans(1:Statistics)=0;
% aggregateWifiTrans(1:Statistics)=0;
% aggregateBmacTrans(1:Statistics)=0;

normalize_w(1:Statistics) = 0;
normalize_b(1:Statistics) = 0;
% tau_w(1:Statistics) = 0;
% tau_b(1:Statistics) = 0;
% Ptr_w(1:Statistics) = 0;
% Ptr_b(1:Statistics) = 0;
xaxis(1:Statistics)=0;
alpha(1:Statistics) = 0;
alpha2(1:Statistics) = 0;
beta(1:Statistics) = 0;
phi(1:Statistics) = 0;

tau(1:Statistics) = 0;
tau2(1:Statistics) = 0;
Pc(1:Statistics) = 0;
aggrtau(1:Statistics) = 0;
aggrPc(1:Statistics) = 0;

runNo=1;

zw_time(1:6)=0;

%TEST TIME SLOTS
% T=100000; 


%Change the following four parameters ONLY!
NODE_WIFI=0;
NODE_BMAC=1;

% WIFI_START=0;
% WIFI_END=0;
% 
% BMAC_START=1;
% BMAC_END=1;

%STATE MACHINE
STATE_CHOOSE=0; % WIFI & BMAC
STATE_INIT=1;
STATE_DIFS=2;
STATE_CSMA=3;
STATE_TRDELAY=4;
STATE_SOFTWARE=5; % for software delay
STATE_OTHER=6; % just make the bmac state not STATE_CHOOSE


%STATE_RETR=4;
%disp(delay(1:N));
cnt=0;
rd=0;
% bmaccount=0;
% wificount=0;
% 
zw_time=clock();
outname=['out-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
    num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
    '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

diary(outname);
diary on;
disp(['*******test start******* at ', ' ', num2str(zw_time(4)),':', ...
    num2str(zw_time(5)),':', num2str(zw_time(6)),' ',num2str(zw_time(2)), ...
    '/', num2str(zw_time(3)),'/', num2str(zw_time(1))]);
diary off;

% tttt='waiting for keyboard input....';

%WIFI
% aMinBE=4;
% aMaxBE=10;
LACK=2;
LACKTIMEOUT=2;
DIFS=6;
OSD_w=30;
% Packet_w=1500;
% aMaxBE = aMinBE*2^5;



%BMAC
slotratio = 3; %slot of boxmac by of wifi
% Wi_b=320;
% Wc_b=80;
OSD_b=220*slotratio;%220*slotratio;
% Packet_b=128;
DELAY_INTERVAL = 1*slotratio;


LDATApb=ceil(Packet_b*8*4/30*slotratio);
LDATAp=ceil(Packet_w*8/54/10);

% for Packet_w=1500:500:1500%[500,1000,1500]
Packet_w=Packet_w+846;
Packet_b=Packet_b+20;
% for aMinBE=6:6
% for aMaxBE=10:10;
% 
% for Packet_b=128:-20:128%[128,108,88,68,48]
% for Wc_b=70:-20:70%[70,50,30]
% for Wi_b=310:-80:310%[310,230,150,70]
LDATAb=ceil(Packet_b*8*4/30*slotratio);
LDATA=ceil(Packet_w*8/54/10);
%for LDATAb=[70,50,30]


for x=1:1 %round
% for Trials=16:16
    diary on;
    zw_time=clock();
    disp(['round: ', num2str(x), ' ', num2str(zw_time(4)),':', num2str(zw_time(5)), ...
        ':', num2str(zw_time(6)),' ',num2str(zw_time(2)),'/', num2str(zw_time(3)), ...
        '/', num2str(zw_time(1))]);
    diary off;
    for wificount=WIFI_START:5:WIFI_END
        for bmaccount=BMAC_START:5:BMAC_END
            for N=bmaccount+wificount:bmaccount+wificount
                if bmaccount==0 && wificount==0
                    break;
                end
                diary on;
                disp('======================================================================');
                disp(['wifi: ', num2str(wificount), ', bmac: ', num2str(bmaccount)]);
                disp(['Wi_b: ', num2str(Wi_b+10), ', Wc_b: ', num2str(Wc_b+10), ', Packet_b: ', num2str(Packet_b), ', LDATA_b: ', num2str(LDATAb), ', OSD_b: ', num2str(OSD_b)]);
                disp(['aMinBE_w: ', num2str(aMinBE), ', aMaxBE_w: ', num2str(aMaxBE), ', Packet_w: ', num2str(Packet_w), ', LDATA_w: ', num2str(LDATA), ', OSD_w: ', num2str(OSD_w), ', DIFS: ', num2str(DIFS)]);
                disp('');
                diary off;
                %time slot=10us
                wifiTRcount(1:TotalTrials)=0;
                bmacTRcount(1:TotalTrials)=0;
                wifiCollision(1:TotalTrials)=0;
                bmacCollision(1:TotalTrials)=0;

                CW(1:TotalTrials)=aMinBE;% aMinBE < BE < aMaxBE
                delay(1:TotalTrials)=0;
                %int32(rand.*(2.^(BE(1:N))));
                %for i=1:N
                %    delay(i)=int32(rand*2^(BE(i)));
                %end

                busyFor(1:TotalTrials)=0; %for transmission delay
                tempbusyFor(1:TotalTrials)=0; %for transmission delay
                nbTransmission(1:TotalTrials)=0;
                nbCollision(1:TotalTrials)=0;
                %nbFailure(1:N)=0;

                DIFScnt(1:TotalTrials)=0; %for DIFS delay
                oneshot(1:TotalTrials)=1; %for oneshot of transmission, 1 means the first try, 0 means the first try fails.
                oneshotcnt(1:TotalTrials)=0; %for oneshot of transmission count of wifi.
                nodeState(1:TotalTrials)=0; %for FSM.
                trfail(1:TotalTrials)=0; %for tr fails.
                nodeType(1:TotalTrials)=0; %for WIFI/BMAC.
                softdelay(1:TotalTrials)=OSD_w; %for WIFI OS delay

                % aMinBEb=3;
                % aMaxBEb=5;
                % macMaxCSMABackoffs=5;
                %NBb(1:N)=0;% <macMaxCSMABackoffs
                CWb(1:TotalTrials)=2;% 
                %BEb(1:N)=aMinBE;% aMinBE < BE < aMaxBE
                nbCCA(1:TotalTrials)=0;
                nbBusy1(1:TotalTrials)=0;
                nbBusy2(1:TotalTrials)=0;

                softwareDelay(1:TotalTrials)=0; %for BMAC OS delay
%                 alpha(1:TotalTrials)=0;
%                 normalize_w(1:TotalTrials) = 0;
%                 normalize_b(1:TotalTrials) = 0;
                %wificount=1;
            %     if N==1
            %         wificount=1;
            %     elseif N==2
            %         wificount=2;
            %     else
            %         wificount=2;
            %         bmaccount=N-2;
            %     end
diary on;

                for i=1:T
                    if mod(i,50000)==0
                        zw_time=clock();
                        disp([num2str(i),' ', num2str(zw_time(4)),':', num2str(zw_time(5)), ...
                            ':', num2str(zw_time(6)),' ',num2str(zw_time(2)),'/', ...
                            num2str(zw_time(3)),'/', num2str(zw_time(1))]);
                        
%                         diary on;
%                         
%                         averageWifiTrans(N)=sum(wifiTRcount-wifiCollision)/(wificount)*Packet_w*8/(T*10/1000000)/1000000;
%                         aggregateWifiTrans(N)=sum(wifiTRcount-wifiCollision)*Packet_w*8/(T*10/1000000)/1000000;
%                         disp(['wifi average throughput: ', num2str(averageWifiTrans(N)), 'Mbps']);
%                         disp(['wifi total throughput: ', num2str(aggregateWifiTrans(N)), 'Mbps']);
%                         averageBmacTrans(N)=sum(bmacTRcount-bmacCollision)/(bmaccount)*Packet_b*8/(T*10/1000000)/1000;
%                         aggregateBmacTrans(N)=sum(bmacTRcount-bmacCollision)*Packet_b*8/(T*10/1000000)/1000;
%                         disp(['bmac average throughput: ', num2str(averageBmacTrans(N)), 'Kbps']);
%                         disp(['bmac total throughput: ', num2str(aggregateBmacTrans(N)), 'Kbps']);
%                     %     disp(['Pcoll: ', num2str(sum(nbCollision)/T)]);
%                         disp(['alpha: ', num2str((alpha(N))/T)]);
%                         disp(['tau_w: ', num2str(sum(wifiTRcount)/(wificount)/T)]);
%                         tau_w(N) = sum(wifiTRcount)/(wificount)/T;
%                         disp(['Ptr_w: ', num2str(sum(wifiTRcount)*(LDATA+LACK)/(wificount)/T)]);
%                         Ptr_w(N) = sum(wifiTRcount)*(LDATA+LACK)/(wificount)/T;
%                         disp(['tau_b: ', num2str(sum(bmacTRcount)/(bmaccount)/T)]);
%                         tau_b(N) = sum(bmacTRcount)/(bmaccount)/T;
%                         disp(['Ptr_b: ', num2str(sum(bmacTRcount)*LDATAb/(bmaccount)/T)]);
%                         Ptr_b(N) = sum(bmacTRcount)*LDATAb/(bmaccount)/T;
%                         disp(['Pcoll_w: ', num2str(sum(wifiCollision)/sum(wifiTRcount))]);
%                         disp(['Pcoll_b: ', num2str(sum(bmacCollision)/sum(bmacTRcount))]);
%                         disp(['normalize_w: ', num2str(sum(normalize_w)/T)]);
%                         disp(['normalize_b: ', num2str(sum(normalize_b)/T)]);
%                         disp(['*', num2str(aggregateWifiTrans(N)), ' ', num2str(aggregateBmacTrans(N))]);
% 
% 
%                         %    disp(averageBmacTrans(N));
%                         disp('======================================================================');
%                         diary off;

                        
                        
                    end

                    if sum(busyFor(1:N))~=0
                        alpha(runNo) = alpha(runNo) + 1;
                    end
                    
                    if sum(busyFor(bmaccount+1:N))~=0
                        normalize_w(runNo) = normalize_w(runNo) + 1;
                    end
                    
                    if sum(busyFor(1:bmaccount))~=0
                        normalize_b(runNo) = normalize_b(runNo) + 1;
                    end

                    for j=1:N
                        if nodeState(j)==STATE_CHOOSE
                            %p=binornd(1,0.8); % binornd(1,x) x is the probability of bmac behavier.
                            %if p==0 % wifi 0.2
                            if j>bmaccount
                                %disp('wifi');
                                nodeType(j)=NODE_WIFI;
                                nodeState(j)=STATE_INIT;
                            else %bmac 0.8
                                %disp('bmac');
                                nodeType(j)=NODE_BMAC;
                                nodeState(j)=STATE_OTHER;
                                %disp('init delay')
                                delay(j)=int32(10+rand*Wi_b)*slotratio;
                                %disp(delay(j));
                            end
                        end
                        %disp('1');
                        if nodeType(j)==NODE_WIFI
                            if nodeState(j)==STATE_INIT
                                DIFScnt(j)=0;
                                if trfail(j)==0
                                    %disp(['[', num2str(j), '] reset oneshot']);
                                    oneshot(j)=0;
                                else
                                    oneshot(j)=0;
                                end
                                nodeState(j)=STATE_DIFS;
                            elseif nodeState(j)==STATE_SOFTWARE
                                softdelay(j)=softdelay(j)-1;
                                if softdelay(j)==0
                                    nodeState(j)=STATE_INIT;
                                    softdelay(j)=OSD_w;
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
                                        
                                        if DIFScnt(j)==DIFS
                                            DIFScnt(j)=0;
                                            if oneshot(j)==1
                                %                     temp = sprintf(' %d', j);
                                %                     temp = sprintf(' %s %d', temp, j);
                                %                     disp(temp);
                                                zw_disp([num2str(i),'  [WIFI:', num2str(j), '] oneshot wifi is sending']);
                                                cnt=cnt+1;
                                                oneshotcnt(j)=oneshotcnt(j)+1;
                                                if cnt>1
                                                    zw_disp([num2str(i),'  [WIFI:', num2str(j), '] wifi collision happen']);
                                                    nbCollision(j)=nbCollision(j)+1;
                                                    nbTransmission(j)=nbTransmission(j)+1;
                                                    wifiCollision(j)=wifiCollision(j)+1;
                                                    wifiTRcount(j)=wifiTRcount(j)+1;
                                                    tempbusyFor(j)=LDATA+LACKTIMEOUT;
                                                    trfail(j)=1;
                                                elseif cnt==1
                                                    rd=j;
                                                end
                                                nodeState(j)=STATE_TRDELAY;
                                            else
                                                if delay(j)==0 %not from freezing
                                                    if trfail(j)==1
%                                                         CW(j)=min(aMaxBE,CW(j)+1);% aMinBE < BE < aMaxBE
%                                                         delay(j)=int32(rand*2^(CW(j)))+2; % must add 1, 
                                                        CW(j)=min(aMaxBE,CW(j)*2);% aMinBE < BE < aMaxBE
                                                        delay(j)=int32(rand*CW(j))+2; % must add 1, 
                                                        zw_disp([num2str(i), '  [WIFI:', num2str(j), '] cw inc to ', num2str(CW(j)), ' backoff for ', num2str(delay(j))]);
                                                    else
%                                                         CW(j)=aMinBE; % aMinBE < BE < aMaxBE
%                                                         delay(j)=int32(rand*2^(CW(j)))+2; % must add 1,
                                                        CW(j)=aMinBE; % aMinBE < BE < aMaxBE
                                                        delay(j)=int32(rand*CW(j))+2; % must add 1,
                                                        zw_disp([num2str(i), '  [WIFI:', num2str(j), '] cw init ', num2str(CW(j)), ' backoff for ', num2str(delay(j))]);
                                                    end

                                                    if delay(j)==0
                                                        cnt=cnt+1;
                                                        zw_disp([num2str(i),'  [WIFI:', num2str(j), '] special wifi is sending...']);
                                                        if cnt>1
                                                            zw_disp([num2str(i),'  [WIFI:', num2str(j), '] wifi collision happen']);
                                                            nbCollision(j)=nbCollision(j)+1;
                                                            nbTransmission(j)=nbTransmission(j)+1;
                                                            wifiCollision(j)=wifiCollision(j)+1;
                                                            wifiTRcount(j)=wifiTRcount(j)+1;
                                                            tempbusyFor(j)=LDATA+LACKTIMEOUT;
                                                            trfail(j)=1;
                                                        elseif cnt==1
                                                            rd=j;
                                                        end
                                                        nodeState(j)=STATE_TRDELAY;
                                                    else
                                                        if sum(busyFor(1:N))==0 % channel idle
                                                            zw_disp([num2str(i),'  [WIFI:', num2str(j), '] <not freezing> backoff for ', num2str(delay(j))]);
                                                            nodeState(j)=STATE_CSMA;
                                                            %delay(j)=delay(j)-1;
                                                            %disp('channel idle');
                                                        else % channle not idle
                                                            zw_disp([num2str(i),'  [WIFI:', num2str(j), '] <not freezing> re-DIFS']);
                                                            oneshot(j)=0;
                                                            DIFScnt(j)=0;
                                                            nodeState(j)=STATE_DIFS;
                                                        end
                                                    end

                                                else %from freezing
                                                    if sum(busyFor(1:N))==0 % channel idle
                                                        zw_disp([num2str(i),'  [WIFI:', num2str(j), '] <freezing> resume backoff for ', num2str(delay(j))]);
                                                        nodeState(j)=STATE_CSMA;
                                                        %delay(j)=delay(j)-1;
                                                        %disp('channel idle');
                                                    else % channle not idle
                                                        zw_disp([num2str(i),'  [WIFI:', num2str(j), '] <freezing> re-DIFS']);
                                                        oneshot(j)=0;
                                                        DIFScnt(j)=0;
                                                        nodeState(j)=STATE_DIFS;
                                                    end      
                                                end
                                            end
                                        end                                        
                                        
                                    else % channle not idle
                                        %disp(['[', num2str(j), '] DIFS un-idle']);
                                        oneshot(j)=0;
                                        DIFScnt(j)=0;
                                    end
                                case{STATE_CSMA}
                                    if sum(busyFor(1:N))==0 % channel idle
                                        nodeState(j)=STATE_CSMA;
                                        delay(j)=delay(j)-1;
                                        if delay(j)==0
                                           zw_disp([num2str(i),'  [WIFI:', num2str(j), '] csma wifi is sending...']);
                                           cnt=cnt+1;
                                            if cnt>1
                                                zw_disp([num2str(i),'  [WIFI:', num2str(j), '] wifi collision happen']);
                                                nbCollision(j)=nbCollision(j)+1;
                                                nbTransmission(j)=nbTransmission(j)+1;
                                                wifiCollision(j)=wifiCollision(j)+1;
                                                wifiTRcount(j)=wifiTRcount(j)+1;
                                                tempbusyFor(j)=LDATA+LACKTIMEOUT;
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
                                       zw_disp([num2str(i),'  [WIFI:', num2str(j), '] wifi send done']);
                                       tempbusyFor(j)=0;
                                        %delay(j)=0;
                                        nodeState(j)=STATE_SOFTWARE;
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
                                if delay(j)==0 || delay(j)==-DELAY_INTERVAL
                                    %disp('cca');
                                    CWb(j)=CWb(j)-1;
                                    if CWb(j)~=0
                                        zw_disp([num2str(i),'  [BMAC:', num2str(j), '] bmac sense 1st idle']);
                                    end
                                    if CWb(j)==0
                                        zw_disp([num2str(i),'  [BMAC:', num2str(j), '] bmac sense 2nd idle, and is sending...']);
                                       cnt=cnt+1;
                                        if cnt>1
                                           zw_disp([num2str(i),'  [BMAC:', num2str(j), '] bmac collision happen']);
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
                                    %zw_disp(['[', num2str(j), '] bmac busyFor(j)=', num2str(busyFor(j))]);
                                    busyFor(j)=busyFor(j)-1;
                                    if busyFor(j)==0

                                        zw_disp([num2str(i),'  [BMAC:', num2str(j), '] bmac send done']);
                                        %NBb(j)=0;% <macMaxCSMABackoffs
                                        CWb(j)=2;% 
                                        %BEb(j)=aMinBEb;% aMinBE < BE < aMaxBE
                                        delay(j)=int32(10+rand*Wi_b)*slotratio+1; % must add 1, 
                                        nodeState(j)=STATE_CHOOSE;
                                        %because this time slot has already used to do busyFor
                                        softwareDelay(j)=OSD_b; % 6.5ms
                                        tempbusyFor(j)=0; %fix bug 6/21/2012
                                    end
                                end
                            end
                       end

                        for j=1:N
                            if nodeType(j)==NODE_BMAC
                                if softwareDelay(j)~=0
                                    continue;
                                end
                                if delay(j)==0 || delay(j)==-DELAY_INTERVAL
                                    if delay(j)==0
                                        nbBusy1(j)=nbBusy1(j)+1;
                                    end
                                    if delay(j)==-DELAY_INTERVAL
                                        nbBusy2(j)=nbBusy2(j)+1;
                                    end                    %disp('backoff');
                    %                 NB(j)=NB(j)+1;% <macMaxCSMABackoffs
                                    CWb(j)=2;% 
                    %                 BE(j)=min(aMaxBE,BE(j)+1);% aMinBE < BE < aMaxBE
                                    %BEb(j)=aMaxBE;% aMinBE < BE < aMaxBE
                                    delay(j)=int32(10+rand*Wc_b)*slotratio+1; % must add 1, 
                                    zw_disp([num2str(i),'  [BMAC:', num2str(j), '] bmac sense busy, backoff for ', num2str(delay(j))]);
                                    %because this time slot has already used to do update
            %                         disp('congestion delay')
            %                         disp(delay(j));
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
                            if delay(j)~=-(DELAY_INTERVAL+1)
                                delay(j)=delay(j)-1;
                            end
                        end
                    end

                    if cnt>1
                        nbCollision(rd)=nbCollision(rd)+1;
                        nbTransmission(rd)=nbTransmission(rd)+1;
                        busyFor=tempbusyFor;
                        if nodeType(rd)==NODE_WIFI
                            zw_disp([num2str(i),'  [WIFI:', num2str(rd), '] wifi collision happen']);
                            wifiCollision(rd)=wifiCollision(rd)+1;
                            wifiTRcount(rd)=wifiTRcount(rd)+1;
                            busyFor(rd)=LDATA+LACKTIMEOUT;
                        else
                            zw_disp([num2str(i),'  [BMAC', num2str(rd), '] bmac collision happen']);
                            bmacCollision(rd)=bmacCollision(rd)+1;
                            bmacTRcount(rd)=bmacTRcount(rd)+1;
                            busyFor(rd)=LDATAb;
                        end
                        trfail(rd)=1;
                    elseif cnt==1
                        nbTransmission(rd)=nbTransmission(rd)+1;
                        if nodeType(rd)==NODE_WIFI
                            wifiTRcount(rd)=wifiTRcount(rd)+1;
                            busyFor(rd)=LDATA+LACK;
                        else
                            bmacTRcount(rd)=bmacTRcount(rd)+1;
                            busyFor(rd)=LDATAb;
                        end
                       trfail(rd)=0;
                    end

                    cnt=0;    

                end
                diary on;

                %disp(wificount);disp(bmaccount);
%                 averageWifiTrans(runNo)=sum(wifiTRcount-wifiCollision)/(wificount)*Packet_w*8/(T*10/1000000)/1000000;
%                 aggregateWifiTrans(runNo)=sum(wifiTRcount-wifiCollision)*Packet_w*8/(T*10/1000000)/1000000;
%                 disp(['wifi average throughput: ', num2str(averageWifiTrans(runNo)), 'Mbps']);
%                 disp(['wifi total throughput: ', num2str(aggregateWifiTrans(runNo)), 'Mbps']);
%                 averageBmacTrans(runNo)=sum(bmacTRcount-bmacCollision)/(bmaccount)*Packet_b*8/(T*10*slotratio/1000000)/1000;
%                 aggregateBmacTrans(runNo)=sum(bmacTRcount-bmacCollision)*Packet_b*8/(T*10*slotratio/1000000)/1000;
%                 disp(['bmac average throughput: ', num2str(averageBmacTrans(runNo)), 'Kbps']);
%                 disp(['bmac total throughput: ', num2str(aggregateBmacTrans(runNo)), 'Kbps']);
                tau(runNo)=sum(wifiTRcount)/N/T*(LDATA);
                tau2(runNo)=sum(wifiTRcount)/N/T;
%                 Pc(runNo)=sum(wifiCollision)/N/T;
%                 aggrtau(runNo)=sum(wifiTRcount)/T;
%                 aggrPc(runNo)=sum(wifiCollision)/T;
                disp(['tau: ', num2str(tau(runNo))]);
                disp(['tau2: ', num2str(tau2(runNo))]);
%                 disp(['Pc: ', num2str(Pc(runNo))]);
%                 disp(['aggrtau: ', num2str(aggrtau(runNo))]);
%                 disp(['aggrPc: ', num2str(aggrPc(runNo))]);
                disp(['Pcoll_w: ', num2str(sum(wifiCollision)/sum(wifiTRcount))]);
%                 disp(['total alpha: ', num2str((alpha(runNo))/T)]);
%                 phi(runNo)=sum(nbCCA)/N/T;
                disp(['phi: ', num2str(sum(nbCCA)/N/T)]);
                disp(['alpha: ', num2str((sum(nbBusy1)/N)/(sum(nbCCA)/N))]);
                alpha2(runNo)=(sum(nbBusy1)/N)/(sum(nbCCA)/N);
                disp(['beta: ', num2str((sum(nbBusy2)/N)/(sum(nbCCA-nbBusy1)/N))]);
                beta(runNo)=(sum(nbBusy2)/N)/(sum(nbCCA-nbBusy1)/N);
                disp(['Pcoll_b: ', num2str(sum(bmacCollision)/sum(bmacTRcount))]);
%                 disp(['tau_w: ', num2str(sum(wifiTRcount)/(wificount)/T)]);
%                 tau_w(runNo) = sum(wifiTRcount)/(wificount)/T;
%                 disp(['Ptr_w: ', num2str(sum(wifiTRcount)*(LDATA+LACK)/(wificount)/T)]);
%                 Ptr_w(runNo) = sum(wifiTRcount)*(LDATA+LACK)/(wificount)/T;
%                 disp(['tau_b: ', num2str(sum(bmacTRcount)/(bmaccount)/T)]);
%                 tau_b(runNo) = sum(bmacTRcount)/(bmaccount)/T;
%                 disp(['Ptr_b: ', num2str(sum(bmacTRcount)*LDATAb/(bmaccount)/T)]);
%                 Ptr_b(runNo) = sum(bmacTRcount)*LDATAb/(bmaccount)/T;
%                 disp(['Pcoll_w: ', num2str(sum(wifiCollision)/sum(wifiTRcount))]);
%                 disp(['Pcoll_b: ', num2str(sum(bmacCollision)/sum(bmacTRcount))]);
                disp(['normalize_all_w: ', num2str((normalize_w(runNo))/T)]);
                disp(['normalize_all_b: ', num2str((normalize_b(runNo))/T)]);
                disp(['normalize_trr_w: ', num2str(sum(wifiTRcount)*(LDATA+LACK)/T)]);
                disp(['normalize_trr_b: ', num2str(sum(bmacTRcount)*LDATAb/T)]);
                disp(['normalize_trrs_w: ', num2str((sum(wifiTRcount)-sum(wifiCollision))*(LDATA+LACK)/T)]);
                disp(['normalize_trrs_b: ', num2str((sum(bmacTRcount)-sum(bmacCollision))*LDATAb/T)]);
                disp(['normalize_trp_w: ', num2str(sum(wifiTRcount)*(LDATAp)/T)]);
                disp(['normalize_trp_b: ', num2str(sum(bmacTRcount)*LDATApb/T)]);
                disp(['normalize_trps_w: ', num2str((sum(wifiTRcount)-sum(wifiCollision))*(LDATAp+LACK)/T)]);
                disp(['normalize_trps_b: ', num2str((sum(bmacTRcount)-sum(bmacCollision))*LDATApb/T)]);
%                 disp(['*', num2str(aggregateWifiTrans(runNo)), ' ', num2str(aggregateBmacTrans(runNo))]);


                %    disp(averageBmacTrans(N));
                disp('======================================================================');
                diary off;
            end
        % end
        % input(tttt);
            xaxis(runNo)=N;
            runNo = runNo+1;

        end
    end
end
% end
% end
% end
% end
% end
% end
%end


%=========for BMAC============
%subplot(221);plot(xaxis, tau,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\tau');axis([0,N,0,1]);%figure(gcf);
%subplot(222);plot(xaxis, Pc,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pc');axis([0,N,0,1]);%figure(gcf);
%subplot(223);plot(xaxis, aggrtau,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('aggr \phi');axis([0,N,0,1]);%figure(gcf);
%subplot(224);plot(xaxis, aggrPc,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('aggr Pc');axis([0,N,0,1]);%figure(gcf);



%=========for BMAC============
% subplot(221);plot(xaxis, phi,'DisplayName','averageCCA','YDataSource','averageCCA');xlabel('Number of nodes');ylabel('\phi');axis([0,N,0,0.03]);%figure(gcf);
% subplot(223);plot(xaxis, alpha2,'DisplayName','averageBusy1','YDataSource','averageBusy1');xlabel('Number of nodes');ylabel('\alpha');axis([0,N,0,1]);%figure(gcf);
% subplot(222);plot(xaxis, beta,'DisplayName','averageBusy2','YDataSource','averageBusy2');xlabel('Number of nodes');ylabel('\beta');axis([0,N,0,1]);%figure(gcf);
% subplot(224);plot(xaxis, alpha,'DisplayName','averageBusy1','YDataSource','averageBusy1');xlabel('Number of nodes');ylabel('\alpha');axis([0,N,0,1]);%figure(gcf);

diary on;
zw_time=clock();
disp(['*******test end******* at ', ' ', num2str(zw_time(4)),':', ...
    num2str(zw_time(5)),':', num2str(zw_time(6)),' ',num2str(zw_time(2)), ...
    '/', num2str(zw_time(3)),'/', num2str(zw_time(1))]);
%     plot(aggregateWifiTrans,'DisplayName','aggregateWifiTrans','YDataSource','aggregateWifiTrans');hold all;
%     plot(aggregateBmacTrans/4,'DisplayName','aggregateBmacTrans','YDataSource','aggregateBmacTrans');hold off;
%     grid off;xlabel('Number of nodes');ylabel('Throughput(Mbps)');axis([0,50,0,100]);
%     figure(gcf);
diary off;

%pause;

ret = 0;
