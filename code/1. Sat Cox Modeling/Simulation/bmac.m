%State Vectors
%N=4;
%time slot=30us

clear all;

Trials=6;
averageCCA(1:Trials)=0;
averageBusy1(1:Trials)=0;
averageBusy2(1:Trials)=0;
averageTrans(1:Trials)=0;
averageSucTrans(1:Trials)=0;

for N=1:Trials


    % aMinBE=3;
    % aMaxBE=5;
    % macMaxCSMABackoffs=5;
    NB(1:Trials)=0;% <macMaxCSMABackoffs
    CW(1:Trials)=2;% 
    % BE(1:N)=aMinBE;% aMinBE < BE < aMaxBE
    cnt=0;
    rd=0;
    length=128; %4ms
    delay(1:Trials)=0;
    %int32(rand.*(2.^(BE(1:N))));
    for i=1:Trials
        delay(i)=int32(10+rand*310); % 310us~10ms
    end
    busyFor(1:Trials)=0;
    T=100000;
    nbTransmission(1:Trials)=0;
    nbCollision(1:Trials)=0;
    nbFailure(1:Trials)=0;
    nbCCA(1:Trials)=0;

    % nodeState(1:N)=0;
    % STATE_INIT=0;
    % STATE_CSMA=1;
    % STATE_TRDELAY=2;
    softwareDelay(1:Trials)=0;


    %disp(delay(1:Trials));

    for i=1:T
        for j=1:N
            if softwareDelay(j)~=0
                continue;
            end
            if delay(j)==0% && lock(j)==0
                nbCCA(j)=nbCCA(j)+1;
            end
        end

        if sum(busyFor(1:N))==0 % channel idle
            %disp('channel idle');
            for j=1:N
                if softwareDelay(j)~=0
                    continue;
                end
                if delay(j)==0 || delay(j)==-7
                    %disp('cca');
                    CW(j)=CW(j)-1;
                    if CW(j)==0
                        cnt=cnt+1;
                        if cnt>1
                            nbCollision(j)=nbCollision(j)+1;
                            nbTransmission(j)=nbTransmission(j)+1;
                            busyFor(j)=length;
                        elseif cnt==1
                            rd=j;
                        end
                    end
                end
            end
            if cnt>1
                nbCollision(rd)=nbCollision(rd)+1;
                nbTransmission(rd)=nbTransmission(rd)+1;
                busyFor(rd)=length;
            elseif cnt==1
                nbTransmission(rd)=nbTransmission(rd)+1;
                busyFor(rd)=length;
            end
        else % channle not idle
            %disp('channel not idle');
            for j=1:N
                if softwareDelay(j)~=0
                    continue;
                end
                if busyFor(j)>0
                    busyFor(j)=busyFor(j)-1;
                    if busyFor(j)==0
                        %disp('sendone');
                        NB(j)=0;% <macMaxCSMABackoffs
                        CW(j)=2;% 
                        %BE(j)=aMinBE;% aMinBE < BE < aMaxBE
                        delay(j)=int32(10+rand*310)+1; % must add 1, 
                        %because this time slot has already used to do busyFor
                        softwareDelay(j)=210; % 6.5ms
                    end
                end
           end

            for j=1:N
                if softwareDelay(j)~=0
                    continue;
                end
                if delay(j)==0 || delay(j)==-7
                    %disp('backoff');
    %                 NB(j)=NB(j)+1;% <macMaxCSMABackoffs
                    CW(j)=2;% 
    %                 BE(j)=min(aMaxBE,BE(j)+1);% aMinBE < BE < aMaxBE
                    %BE(j)=aMaxBE;% aMinBE < BE < aMaxBE
                    delay(j)=int32(10+rand*80)+1; % must add 1, 
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
        cnt=0;

        for j=1:N
            if softwareDelay(j)~=0
                softwareDelay(j)=softwareDelay(j)-1;
                continue;
            end
            if delay(j)~=-8
                delay(j)=delay(j)-1;
            end
        end

    end
    
    
    averageTrans(N)=sum(nbTransmission-nbCollision)/N*8*100/(T*30/1000000)/1000;
    disp(averageTrans(N));

    
    
end
plot(averageTrans,'DisplayName','averageTrans','YDataSource','averageTrans');xlabel('Number of nodes');ylabel('Throughput(Kbps)');axis([0,50,0,60]);figure(gcf);
% for i=1:T
%     for j=1:N
%         if nodeState(j)==STATE_INIT
%             delay(j)=int32(rand*2^(BE(j)));
%             nodeState(j)=STATE_CSMA;
%         end
%         
%         switch nodeState(j)
%             case{STATE_CSMA}
%                 
%         end
%         
%         if delay(j)==0% && lock(j)==0
%             nbCCA(j)=nbCCA(j)+1;
%         end
%         
%     end
% end

