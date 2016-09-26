% m=4;
% n=3;
clear all;
format shortG;

runNo=1;

WIFI_START=10;
WIFI_END=10;

BMAC_START=50;
BMAC_END=50;
STEP_SIZE=5;

totalRuns=(WIFI_END-WIFI_START)/STEP_SIZE+1;

xaxis(1:totalRuns)=0;
averThr_w(1:totalRuns)=0;
aggrThr_w(1:totalRuns)=0;
averThr_b(1:totalRuns)=0;
aggrThr_b(1:totalRuns)=0;
alphas(1:totalRuns)=0;
tau_ws(1:totalRuns)=0;
tau_bs(1:totalRuns)=0;

Wis(1:totalRuns)=0;
Wcs(1:totalRuns)=0;
Eplbs(1:totalRuns)=0;
Wis(1:totalRuns)=0;
Wcs(1:totalRuns)=0;
Eplbs(1:totalRuns)=0;
m0s(1:totalRuns)=0;
ms(1:totalRuns)=0;
Eplws(1:totalRuns)=0;

taus(1:totalRuns)=0;
Pcs(1:totalRuns)=0;


global N;
% global Ls;
global Lo;
global CW;
global m;
global bolimit;
global packetsize;
global windowsize;
global nn;

global L_s;
global L_c;
global L_sc;

global W_d;
global m_d;
global l_d;
global s_d;
global w_d;

global W_a;
global m_a;
global l_a;
global s_a;
global w_a;

disp('**************************Markov Chain Model***************************');
disp('**************************parameters***************************');

%for ATIM window
%=========================================
Lo=0;%20;

% W=32;
% m=6; %m is the CW max.
% l=9; %j is the retry limit.
% s = 6; %packet size (slots)
% a = 100; %ATIM window size (slots)


% W_a=32; %contention window size
% m_a=5; %the CW max.
% l_a=7; %the retry limit.
% s_a = 3; %ATIM packet size (slots)
% w_a = 100; %ATIM window size (slots)

W_a=3; %contention window size
m_a=0; %the CW max.
l_a=0; %the retry limit.
s_a = 2; %ATIM packet size (slots)
w_a = 5; %ATIM window size (slots)

nn = 3; %number of total nodes, including myself
% no = 4; %number of active nodes, including myself (only) if I am active
%=========================================

% test = 0;
% for i=1:l+1
%     if i<m+1
%         test = test + W*2^(i-1);%cw length
%      else %i>=m+1
%         test = test + W*2^m;%cw length
%     end
% end

% if test>a
%     error(['test=' num2str(test) ' must be greater than a!']);
%     return;
% end


w0 = CW; %contention window size (slots)
fprintf('# nodes: %d\n', nn);
fprintf('====ATIM window related:====\n');
fprintf('CW: %d, max: %d, limit: %d\n', W_a, m_a, l_a);
fprintf('packet size: %d, atim: %d\n', s_a, w_a);







%for DATA window
%=========================================
w_d=100;%20; %DATA window size (slots)
s_d=10;%4; %DATA packet size (slots)
W_d=W_a;% con
m_d=m_a; %m is the CW max.
l_d=l_a; %j is the retry limit.

%=========================================

L_s=36;
L_c=38;
L_sc=37;
Payload = L_s - 16;

fprintf('====DATA window related:====\n');
fprintf('CW: %d, max: %d, limit: %d\n', W_d, m_d, l_d);
fprintf('packet size: %d, atim: %d\n', s_d, w_d);

disp('******************************ATIM start************************************');
SYMBOL_VERIFICATION = 0;
global snarray;
global cwarray;

for i=1:1
    
%     T_matrix=[1 2 7 0 0; 2 4 4 0 0; 3 5 2 0 0 ; 0 0 0 0 10; 0 0 10 0 0];
%     T_matrix=T_matrix/10;
%     [K, ~]=size(T_matrix);
%     
%     Q_matrix=sparse(1:K,1:K,1)-T_matrix.';
%     
%     disp(Q_matrix);
%     
%     b=[zeros(K-1,1);1];
%     
%     Q_matrix(K,:)=1;
%     disp(Q_matrix);
%     
%     pi=Q_matrix\b;
%     
%     disp(pi);
%     
%     return;
    
    windowsize=w_a;%20; %ATIM window size (slots)
    packetsize=s_a;%4; %ATIM packet size (slots)
    CW=W_a;
    bolimit=l_a;
    m=m_a;
    count=0;
    
    snarray=[];
    cwarray=[];
    cwsave=0;
    for a=1:bolimit+1
        if a<m+1
            cw = CW*2^(a-1);
        else
            cw = CW*2^m;
        end
        cwarray=[cwarray cw];
%         disp(cw);
%         for b=cw:-1:1
%             for c=windowsize:-1:1
%                 for d=nn:-1:1
% %                     fprintf('a,b,c,d=%d,%d,%d,%d\n', a, b, c, d);
%                     count=count+1;
%                 end
%             end
%         end
        count = count+cw*windowsize*nn;
        snarray=[snarray count];
%         cwsave=cwsave+count;
        
%         fprintf('count:%d\n', count);
    end
    
    disp(cwarray);
    disp(snarray);
    
    
    
    test=1;
    
%     no=activenn(test);
    for test=31:1:31
        fprintf('test: %d ,[%d,%d,%d,%d]\n', test, bostage(test)-1, bocounter(test)-1, windowleft(test), activenn(test));
        temp = targets(test);
        test1 = temp(1);
        fprintf('\t\t first: %d, [%d,%d,%d,%d]\n', test1, bostage(test1)-1, bocounter(test1)-1, windowleft(test1), activenn(test1));
        test1 = temp(2);
        fprintf('\t\t second: %d, [%d,%d,%d,%d]\n', test1, bostage(test1)-1, bocounter(test1)-1, windowleft(test1), activenn(test1));
        test1 = temp(3);
        fprintf('\t\t third: %d, [%d,%d,%d,%d]\n\n', test1, bostage(test1)-1, bocounter(test1)-1, windowleft(test1), activenn(test1));
        
%         fprintf('ativenn: %d ,', activenn(test));
%         fprintf('windowleft: %d ,', windowleft(test));
%         fprintf('bocounter: %d ,', bocounter(test));
%         fprintf('bostage: %d \n', bostage(test));
%         targets(test)
    end
    
%     return;
    
    matrix_size = 1;%count*(1+2*packetsize) + nn*windowsize;
    
    T_matrix = sparse(1,1,0,matrix_size,matrix_size);
    
    tau=0.4;
    tauc=1-tau;

    p_a=@(no) tauc^(no-1);
    p_c=@(no) (no-1)*tau*tauc^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
    p_b=@(no) 1-p_a(no)-p_c(no);
    p_d=@(no) p_b(no)+p_c(no);    
    
    
    windowendsstate = count+1;
    
    waitwinends_offset = count+1;
    
    fprintf('waitwinends_offset: %d\n', waitwinends_offset);
    
    txstate_offset=count+1+(windowsize-packetsize)*nn;
    for a=1:count
%         disp(a);
        temp = targets(a);
        disp(temp);
        if temp(1)~=-1 && temp(1)~=0 && temp(1)~=-2
            T_matrix(a,temp(1))=p_a(activenn(a));
        elseif temp(1)==-1
            
            T_matrix(a,windowendsstate)=1;
            continue;
            
%             T_matrix(a,offset)=p_a(activenn(a));
% %             disp(offset);
%             offset=offset+1;
% %             fprintf('windowleft: %d\n', windowleft(a));
%             for b=2:windowleft(a)-1
%                 T_matrix(offset-1,offset)=1;
% %                 disp(offset);
%                 offset=offset+1;
%                 %T_matrix(a,temp(2))=p_b(activenn(a));
%             end
% %             offset=offset+1;
%             T_matrix(offset-1,windowendsstate)=1;
% %             disp(offset);
        elseif temp(1)==-2% i.e. bocounter(a)==1
            T_matrix(a,waitwinends_offset+mod(a,windowsize*nn))=p_a(activenn(a));
        end
        if temp(2)~=-1 && temp(2)~=0 && temp(2)~=-2
            T_matrix(a,txstate_offset)=p_b(activenn(a));
%             disp(offset);
            txstate_offset=txstate_offset+1;
            for b=2:packetsize-1
                T_matrix(txstate_offset-1,txstate_offset)=1;
%                 disp(offset);
                txstate_offset=txstate_offset+1;
                %T_matrix(a,temp(2))=p_b(activenn(a));
            end
%             offset=offset+1;
            T_matrix(txstate_offset-1,temp(2))=1;
%                 disp(offset);
        elseif temp(2)==-1
            T_matrix(a,txstate_offset)=p_b(activenn(a));
%             disp(offset);
            txstate_offset=txstate_offset+1;
%             fprintf('windowleft: %d\n', windowleft(a));
            for b=2:windowleft(a)-1
                T_matrix(txstate_offset-1,txstate_offset)=1;
%                 disp(offset);
                txstate_offset=txstate_offset+1;
                %T_matrix(a,temp(2))=p_b(activenn(a));
            end
%             offset=offset+1;
            T_matrix(txstate_offset-1,windowendsstate)=1;
%             disp(offset);
        elseif temp(2)==-2
        end
        if temp(3)~=-1 && temp(3)~=0 && temp(3)~=-2
            T_matrix(a,txstate_offset)=p_c(activenn(a));
%             disp(offset);
            txstate_offset=txstate_offset+1;
            for b=2:packetsize-1
                T_matrix(txstate_offset-1,txstate_offset)=1;
%                 disp(offset);
                txstate_offset=txstate_offset+1;
                %T_matrix(a,temp(2))=p_b(activenn(a));
            end
%             offset=offset+1;
            T_matrix(txstate_offset-1,temp(3))=1;
%                 disp(offset);
%             T_matrix(a,temp(3))=p_c(activenn(a));
        elseif temp(3)==-1
            T_matrix(a,txstate_offset)=p_c(activenn(a));
%             disp(offset);
            txstate_offset=txstate_offset+1;
            for b=2:windowleft(a)-1
                T_matrix(txstate_offset-1,txstate_offset)=1;
%                 disp(offset);
                txstate_offset=txstate_offset+1;
                %T_matrix(a,temp(2))=p_b(activenn(a));
            end
%             offset=offset+1;
            T_matrix(txstate_offset-1,windowendsstate)=1;
%             disp(offset);
            
        elseif temp(3)==-2
            
        end
   
    end
%     
%     disp(T_matrix);
    spy(T_matrix);
    
    
%     for a=1:size(T_matrix,1)
%         
%         if sum(T_matrix(a,:))~=1 && a~=windowendsstate
%             disp(a);
%             error('error!');
%         end
%         
%     end
    
    for a=1:size(T_matrix,1)
        [~,ss] = size(find(T_matrix(a,:)~=0));
        tt=find(T_matrix(a,:)~=0);
        for b=1:ss
            fprintf('(%d,%d) %g\n', a, tt(b), (full(T_matrix(a,tt(b)))));
%             disp(full(T_matrix(a,tt(b))));
        end
        fprintf('\n');
    end
    
    
    
    
%     test=2;
%     
% 
%     
%     disp(test);
    
    
    
    
    return;
    
    
    
    
    
    
    
    
    
    
    
    [K, ~]=size(T_matrix);
    
    Q_matrix=sparse(1:K,1:K,1)-T_matrix.';
    
    disp(Q_matrix);
    
    b=[zeros(K-1,1);1];
    
    Q_matrix(K,:)=1;
    disp(Q_matrix);
    
    pi=Q_matrix\b;
    
    disp(pi);
    
    %   T_matrix=[1 2 7 ; 2 4 4  ; 3 5 2 ];
    % T_matrix=T_matrix/10;
    %  [K, ~]=size(T_matrix);
    %
    %  Q_matrix=sparse(1:K,1:K,1)-T_matrix.';
    %
    %  disp(Q_matrix);
    %
    % b=[zeros(K-1,1);1];
    %
    %     Q_matrix(K,:)=1;
    %     disp(Q_matrix);
    %
    %     pi=Q_matrix\b;
    %
    %     disp(pi);
    %
    %
    % T_matrix=[1 2 7 0 ; 2 4 0 4 ; 3 5 2 0  ; 0 0 10 0 ];
    % T_matrix=T_matrix/10;
    %  [K, ~]=size(T_matrix);
    %
    %  Q_matrix=sparse(1:K,1:K,1)-T_matrix.';
    %
    %  disp(Q_matrix);
    %
    % b=[zeros(K-1,1);1];
    %
    %     Q_matrix(K,:)=1;
    %     disp(Q_matrix);
    %
    %     pi=Q_matrix\b;
    %
    %     disp(pi);
    
    return;
    
    
    xguess=[0.2 0.2 0.1]';
    options = optimset('MaxFunEvals',2000, 'MaxIter', 2000);
    xvect = fsolve('solver4d_unsat_math', xguess, options);
    tau = xvect(1);
    tau2 = xvect(2);
    init = xvect(3);
    
    %     tau = xvect(3);
    
    disp('======================================================================');
    disp(['  N = ', num2str(N) ])
    
    disp('The roots from the default "fsolve" are: ')
    %     disp(['  p = ', num2str(p) ]);
    fprintf('  init: %g\n', init);
    
    %         zw0 = 0;
    %         zw1 = 0;
    %         for i=1:l+1
    %
    %             zw0 = zw0 + (sum(atim_backoff_states_container{1,i}(:,1)));
    %             zw1 = zw1 + (sum(sum(atim_backoff_states_container{1,i})));
    %
    %
    %         end
    %             tau = zw0/zw1;
    
    disp(['  tau = ', num2str(tau) ]);
    disp(['  tau2 = ', num2str(tau2) ]);
    
    taus(runNo)=tau;
    %     Pcs(runNo)=p;
    
    aggrThr_w(runNo)=N*tau*(1-tau)^(N-1)*Payload/((1-tau)^N+N*tau*(1-tau)^(N-1)*L_s+(1-(1-tau)^N-N*tau*(1-tau)^(N-1))*L_c);
    disp(['  throughput = ', num2str(aggrThr_w(runNo)) ]);
    
    
    xaxis(runNo)=N;
    
    runNo = runNo+1;
end

subplot(221);plot(xaxis, taus,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\tau');axis([0,N,0,0.1]);hold all;%figure(gcf);
subplot(222);plot(xaxis, Pcs,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pc');axis([0,N,0,1]);hold all;%figure(gcf);
subplot(223);plot(xaxis, aggrThr_w,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('throughput');axis([0,N,0,1]);hold all;%figure(gcf);



