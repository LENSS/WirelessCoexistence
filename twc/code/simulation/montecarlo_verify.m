%State Vectors
%N=5;

%disp(delay(1:N));

Trials=50;
averageCCA(1:Trials)=0;
averageBusy1(1:Trials)=0;
averageBusy2(1:Trials)=0;
averageTrans(1:Trials)=0;
averageSucTrans(1:Trials)=0;
T=100000;
aMinBE=3;
aMaxBE=5;
macMaxCSMABackoffs=5;
cnt=0;
rd=0;
length=15;

for N=1:Trials
    
NB(1:Trials)=0;% <macMaxCSMABackoffs
CW(1:Trials)=2;% 
BE(1:Trials)=aMinBE;% aMinBE < BE < aMaxBE
%delay(1:N)=1;
%int32(rand.*(2.^(BE(1:N))));
busyFor(1:Trials)=0;
nbTransmission(1:Trials)=0;
nbCollision(1:Trials)=0;
nbFailure(1:Trials)=0;
nbCCA(1:Trials)=0;
nbBusy1(1:Trials)=0;
nbBusy2(1:Trials)=0;

    for i=1:N
        delay(i)=int32(rand*2^(BE(i)));
    end

    for i=1:T
        for j=1:N
            if delay(j)==0% && lock(j)==0
                nbCCA(j)=nbCCA(j)+1;
            end
        end

        if sum(busyFor(1:N))==0 % channel idle
            %disp('channel idle');
            for j=1:N
                if delay(j)==0 || delay(j)==-1
                    %disp('cca');
                    CW(j)=CW(j)-1;
                    if CW(j)==0
                        cnt=cnt+1;
                        if cnt>1
                            nbCollision(j)=nbCollision(j)+1;
                            nbTransmission(j)=nbTransmission(j)+1;
                            busyFor(j)=busyFor(j)+length;
                        elseif cnt==1
                            rd=j;
                        end
                    end
                end
            end
            if cnt>1
                nbCollision(rd)=nbCollision(rd)+1;
                nbTransmission(rd)=nbTransmission(rd)+1;
                busyFor(rd)=busyFor(rd)+length;
            elseif cnt==1
                nbTransmission(rd)=nbTransmission(rd)+1;
                busyFor(rd)=busyFor(rd)+length;
            end
        else % channle not idle
            %disp('channel not idle');
            for j=1:N
                if busyFor(j)>0
                    busyFor(j)=busyFor(j)-1;
                    if busyFor(j)==0
                        %disp('sendone');
                        NB(j)=0;% <macMaxCSMABackoffs
                        CW(j)=2;% 
                        BE(j)=aMinBE;% aMinBE < BE < aMaxBE
                        delay(j)=int32(rand*2^(BE(j)))+1; % must add 1, 
                        %because this time slot has already used to do busyFor
                    end
                end
           end

            for j=1:N
                if delay(j)==0 || delay(j)==-1
                    if delay(j)==0
                        nbBusy1(j)=nbBusy1(j)+1;
                    end
                    if delay(j)==-1
                        nbBusy2(j)=nbBusy2(j)+1;
                    end                    %disp('backoff');
                    NB(j)=NB(j)+1;% <macMaxCSMABackoffs
                    CW(j)=2;% 
                    BE(j)=min(aMaxBE,BE(j)+1);% aMinBE < BE < aMaxBE
                    delay(j)=int32(rand*2^(BE(j)))+1; % must add 1, 
                    %because this time slot has already used to do update
                    if NB(j)==macMaxCSMABackoffs
                        nbFailure(j)=nbFailure(j)+1;
                        NB(j)=0;% <macMaxCSMABackoffs
                        CW(j)=2;% 
                        BE(j)=aMinBE;% aMinBE < BE < aMaxBE
                        delay(j)=int32(rand*2^(BE(j)))+1; % must add 1, 
                        %because this time slot has already used to do reset
                    end
                end
            end

        end
        cnt=0;

        for j=1:N
            if delay(j)~=-2
                delay(j)=delay(j)-1;
            end
        end

    end
    %disp(nbCCA)
    disp(sum(nbCCA)/N/T);
    averageCCA(N)=sum(nbCCA)/N/T;
    averageBusy1(N)=(sum(nbBusy1)/N)/(sum(nbCCA)/N);
    averageBusy2(N)=(sum(nbBusy2)/N)/(sum(nbCCA-nbBusy1)/N);
    averageTrans(N)=(sum(nbTransmission)/N/T);
    averageSucTrans(N)=(sum(nbTransmission-nbCollision)/N/T);

    %averageCCA(N)=sum(nbCCA)/N/T;
    
end


subplot(221);plot(averageCCA,'DisplayName','averageCCA','YDataSource','averageCCA');xlabel('Number of nodes');ylabel('\phi');axis([0,50,0,0.09]);figure(gcf);
subplot(223);plot(averageBusy1,'DisplayName','averageBusy1','YDataSource','averageBusy1');xlabel('Number of nodes');ylabel('\alpha');axis([0,50,0,1]);figure(gcf);
subplot(222);plot(averageBusy2,'DisplayName','averageBusy2','YDataSource','averageBusy2');xlabel('Number of nodes');ylabel('\beta');axis([0,50,0,0.9]);figure(gcf);
subplot(224);plot(averageSucTrans,'DisplayName','averageSucTrans','YDataSource','averageSucTrans');hold all;plot(averageTrans,'DisplayName','averageTrans','YDataSource','averageTrans');xlabel('Number of nodes');ylabel('P_s');axis([0,50,0,0.02]);hold off;figure(gcf);


