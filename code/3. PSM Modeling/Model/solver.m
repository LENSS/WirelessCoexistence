% m=4;
% n=3;
clear all;

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
global W;
global m;
global l;
global s;
global a;
global nn;

global L_s;
global L_c;
global L_sc;


disp('**************************Markov Chain Model***************************');
disp('******************************start************************************');

%for ATIM window
%=========================================
Lo=0;%20;

W=32;
m=6; %m is the CW max.
l=9; %j is the retry limit.
s = 6; %packet size (slots)
a = 100; %ATIM window size (slots)

% W=12;
% m=6; %m is the CW max.
% l=6; %j is the retry limit.
% s = 2; %packet size (slots)
% a = 32; %ATIM window size (slots)

nn = 1; %meanningless here
no = 4; %nodes number, including myself
%=========================================

w0 = W; %contention window size (slots)


%for DATA window
%=========================================
d = 20;
t = 4;

%=========================================

L_s=36;
L_c=38;
L_sc=37;
Payload = L_s - 16;

global atim_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global atim_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_atim_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_atim_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container

for i=1:l+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global A',num2str(i),';']);
    eval(['global B',num2str(i),';']);
%     eval(['A',num2str(i),'=zeros(1,W*2^(i-1))']);
%     B{i}=eval(['A',num2str(i)]);
end

global wait_atimwin_ends; %waiting for the end of the ATIM window
global wait_atimwin_ends_temp;
global sum_wait_atimwin_ends; %for storing the sum of wait_atimwin_ends

% wait_atimwin_ends = sym(zeros(a-s+1,1)); %store the ATIM waiting chain
% wait_atimwin_ends_temp = sym(zeros(a-s+1,1)); %store the ATIM waiting chain temporarily
wait_atimwin_ends = (zeros(a-s+1,1)); %store the ATIM waiting chain
wait_atimwin_ends_temp = (zeros(a-s+1,1)); %store the ATIM waiting chain temporarily

global wait_datawin_ends;
global wait_datawin_ends_temp;
global sum_wait_datawin_ends; %for storing the sum of wait_datawin_ends

wait_datawin_ends = sym(zeros(d-t+1,1)); %store the DATA waiting chain
wait_datawin_ends_temp = sym(zeros(d-t+1,1)); %store the DATA waiting chain



% return;
    
for N=WIFI_START:5:WIFI_END

%     Ls=0;%40;


%     A = sym(zeros(k,W));
%     A = (zeros(k,W));
%     
%     return;
    
    for i=1:l+1
        if i<m+1
%             eval(['A',num2str(i),'=sym(zeros(a*nn,W*2^(i-1)));']);
%             eval(['B',num2str(i),'=sym(zeros(s*nn,W*2^(i-1)));']);
            eval(['A',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
            eval(['B',num2str(i),'=zeros(s*nn,W*2^(i-1));']);
        else
%              eval(['A',num2str(i),'=sym(zeros(a*nn,W*2^m));']);
%              eval(['B',num2str(i),'=sym(zeros(s*nn,W*2^m));']);
             eval(['A',num2str(i),'=zeros(a*nn,W*2^m);']);
             eval(['B',num2str(i),'=zeros(s*nn,W*2^m);']);
        end
        atim_backoff_states_container{i}=eval(['A',num2str(i)]);
        atim_windowends_states_container{i}=eval(['B',num2str(i)]);
%        disp(atim_backoff_states_container{1,i});
%        disp(atim_windowends_states_container{1,i});
    end
    
%     %     return;
%     
% %     syms p;
%     p=0.2;
%     sum_atim_backoff_states = 0;
%     sum_atim_windowends_states = 0;
%     sum_wait_atimwin_ends = 0;
% 
%     for i=1:l+1
%         
%         if i<m+1
%             k = W*2^(i-1);%cw length
%          else %i>=m+1
%             k = W*2^m;%cw length
%             
%         end
%            
%             if i==1
%                 for j=k:-1:1
%                     atim_backoff_states_container{1,i}(1,j) = 1/k;
%                 end
%             else
%                 Ct=atim_backoff_states_container{1,i-1}(:,1);
%                 
% %                 disp(simplify(Ct));
% 
%                 for j=k:-1:1
%                     for n=1:a*nn
%                         if Ct(n,1)~=0
%                             if n+s*nn<=a*nn
%                                 atim_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p/k+atim_backoff_states_container{1,i}(n+s*nn,j);
%                             end
%                         end
%                     end
%                 end            
%             end
%             
% %             
% %             disp('atim_backoff_states_container{1,i}');
% %             disp(atim_backoff_states_container{1,i});
%             
%            for j=k:-1:1
%                 
%                 for n=1:a*nn
%                 %             disp(m);
%                 %             disp(n);
%                     if atim_backoff_states_container{1,i}(n,j)~=0
%                         
%                         if j~=1
%                             if n<=a*nn-s+1
%                                 sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p);
%                                 %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
%                             else
%                                 sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a*nn-n)*p);
%                                 %here since if ATIM window ends, the freezing time may not be s-1 long.
%                             end
%                         else %j==1
%                             if n<=a*nn-s+1
%                                 sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1));
%                             else
%                                 sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a*nn-n)*1);
%                             end
%                         end
% %                         disp('sum_atim_backoff_states:');
% %                         disp(simplify(sum_atim_backoff_states));
%                         
%                         if j>1 %except the first column
%                 %                     disp('m>1');
%                             if n+1*nn<=a*nn
%                                 atim_backoff_states_container{1,i}(n+1*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*(1-p)+atim_backoff_states_container{1,i}(n+1*nn,j-1);
%                             end          
%                             if n+s*nn<=a*nn
%                                 atim_backoff_states_container{1,i}(n+s*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p+atim_backoff_states_container{1,i}(n+s*nn,j-1);
%                             end
%                         end
% 
%                         if n+1*nn>a*nn
%                             atim_windowends_states_container{1,i}(n+1*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*(1-p)+atim_windowends_states_container{1,i}(n+1*nn-a*nn,j);
%             %                         disp(Bt(m-1,n+1-k));
%                         end
% 
%                         if n+s*nn>a*nn
%                             atim_windowends_states_container{1,i}(n+s*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p+atim_windowends_states_container{1,i}(n+s*nn-a*nn,j);
%             %                         disp(Bt(m-1,n+s-k));
%                         end
%                     end
%                 end
%            end
%            
% %             disp('atim_backoff_states_container{1,i} after...');
% % %             disp(simplify(atim_backoff_states_container{1,i}));
% %             disp((atim_backoff_states_container{1,i}));
% %             
% %             disp('atim_windowends_states_container{1,i} after...');
% % %             disp(simplify(atim_windowends_states_container{1,i}));
% %             disp((atim_windowends_states_container{1,i}));
% 
% %             
% %             disp('sum_atim_backoff_states:');
% % %             sum_atim_backoff_states = simplify(sum_atim_backoff_states);
% %             disp(simplify(sum_atim_backoff_states));
%             
% %             break;
%         
% %         syms a1 b1 c1 d1 e1 f1 g1 h1 d2 e2 f2 g2 h2;
% %         zw = sym([a1;b1;c1;d1;e1;f1;g1;h1]);
% %         zw1 = sym([0;0;0;d2;e2;f2;g2;h2]);
% % %         disp(zw);
% 
% 
% %         wait_atimwin_ends_temp = sym(zeros(a-s+1,1)); %reset wait_atimwin_ends_temp
%         wait_atimwin_ends_temp = zeros(a-s+1,1); %reset wait_atimwin_ends_temp
%         
%         
% % 
% % 
% %          for n=1:8
% %         %             disp(m);
% %         %             disp(n);
% %             if zw(n,1)~=0
% %                 if n==1
% %                     wait_atimwin_ends_temp(n,1) = zw(n,1);
% %                 else
% %                     wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1,1)+ zw(n,1);
% %                 end
% %             end
% %          end
% %          
% %         disp('wait_atimwin_ends_temp:');
% %         disp(simplify(wait_atimwin_ends_temp));
% % 
% %         
% %         wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
% %         disp('wait_atimwin_ends:');
% %         disp(simplify(wait_atimwin_ends));
% %         wait_atimwin_ends_temp = sym(zeros(a-s+1,1)); %reset wait_atimwin_ends_temp
% %         
% %          for n=1:8
% %         %             disp(m);
% %         %             disp(n);
% %             if zw1(n,1)~=0
% %                 if n==1
% %                     wait_atimwin_ends_temp(n,1) = zw1(n,1);
% %                 else
% %                     wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1,1)+ zw1(n,1);
% %                 end
% %             end
% %          end
% %          
% %         disp('wait_atimwin_ends_temp:');
% %         disp(simplify(wait_atimwin_ends_temp));
% % 
% %         
% %         wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
% %         disp('wait_atimwin_ends:');
% %         disp(simplify(wait_atimwin_ends));
% % 
% % return;
% % 
% 
%          for n=1:a*nn-s+1
%         %             disp(m);
%         %             disp(n);
%             if atim_backoff_states_container{1,i}(n,1)~=0
%                 if n==1
%                     wait_atimwin_ends_temp(n,1) = atim_backoff_states_container{1,i}(n,1)*(1-p);
%                 else
%                     wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1,1)+ atim_backoff_states_container{1,i}(n,1)*(1-p);
%                 end
%             end
%          end
%          
% %         disp('wait_atimwin_ends_temp:');
% %         disp(simplify(wait_atimwin_ends_temp));
% 
%         
%         wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
%         
% %         disp('atim_backoff_states_container:');
% %         disp(simplify(atim_backoff_states_container{1,i}));
% %         disp('atim_windowends_states_container:');
% %         disp(simplify(atim_windowends_states_container{1,i}));
% %         
% %         
%         sum_atim_windowends_states = sum_atim_windowends_states + sum(sum(atim_windowends_states_container{1,i}));
% %         disp('sum_atim_windowends_states:');
% %         disp(sum_atim_windowends_states);
% %         
% %         disp('wait_atimwin_ends:');
% %         disp(simplify(wait_atimwin_ends));
%     end
% 
%         disp('sum_atim_backoff_states:');
% %         sum_atim_backoff_states = simplify(sum_atim_backoff_states);
%         disp(sum_atim_backoff_states);   
%         
%         disp('sum_atim_windowends_states:');
% %         sum_atim_windowends_states = simplify(sum_atim_windowends_states);
%         disp(sum_atim_windowends_states);
%         
% %         sum_wait_atimwin_ends = simplify(sum(wait_atimwin_ends));
%         sum_wait_atimwin_ends = sum(wait_atimwin_ends);
%         disp('sum_wait_atimwin_ends:');
%         disp(sum_wait_atimwin_ends);
%         
%         
%         all = sum_atim_backoff_states + sum_atim_windowends_states + sum_wait_atimwin_ends;
%         fprintf('all: %g\n', all);
%         
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
%             fprintf('tau: %g\n', tau);
%     return;
% 

   
%     syms init p;
%     for j=1:m
%         for i=W*2^(j-1):-1:1
%             B{1,j}(i) = B{1,j}(i) + init*p^(j-1)/(W*2^(j-1));
%             if i>1
%                 B{1,j}(i-1) = B{1,j}(i-1) + B{1,j}(i);
%             end
%         end
%         disp(B{1,j});
%     end
%         for i=W*2^m:-1:1
%             B{1,m+1}(i) = B{1,m+1}(i) + B{1,m}(1)*p/(1-p)/(W*2^(m));
%              if i>1
%                 B{1,m+1}(i-1) = B{1,m+1}(i-1) + B{1,m+1}(i);
%              end
%         end
%         disp(B{1,m+1});
%     
%     for j=1:m+1
%         sum(B{1,j})
%     end
%     s = 0;
%     for j=1:m+1
%         s = s + sum(B{1,j});
%     end
%         ss = simplify(s);
%         disp(ss);
%         
%     sss = solve(ss-1,'init');
%     disp(sss/(1-p));
%     
%     return;
    
    
    xguess=[0.8 0.1]';
    options = optimset('MaxFunEvals',2000, 'MaxIter', 2000);
    xvect = fsolve('solver_math', xguess, options);
    p = xvect(1);
    init = xvect(2);
    
%     tau = xvect(3);

    disp('======================================================================');
    disp(['  N = ', num2str(N) ])

    disp('The roots from the default "fsolve" are: ')
    disp(['  p = ', num2str(p) ]);
    fprintf('  init: %g\n', init);
    
        zw0 = 0;
        zw1 = 0;
        for i=1:l+1
            
            zw0 = zw0 + (sum(atim_backoff_states_container{1,i}(:,1)));
            zw1 = zw1 + (sum(sum(atim_backoff_states_container{1,i})));
            
            
        end
            tau = zw0/zw1;
            
    disp(['  tau = ', num2str(tau) ]);

    taus(runNo)=tau;
    Pcs(runNo)=p;
    
    aggrThr_w(runNo)=N*tau*(1-tau)^(N-1)*Payload/((1-tau)^N+N*tau*(1-tau)^(N-1)*L_s+(1-(1-tau)^N-N*tau*(1-tau)^(N-1))*L_c);
    disp(['  throughput = ', num2str(aggrThr_w(runNo)) ]);
        

    xaxis(runNo)=N;

    runNo = runNo+1;
end

subplot(221);plot(xaxis, taus,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\tau');axis([0,N,0,0.1]);hold all;%figure(gcf);
subplot(222);plot(xaxis, Pcs,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pc');axis([0,N,0,1]);hold all;%figure(gcf);
subplot(223);plot(xaxis, aggrThr_w,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('throughput');axis([0,N,0,1]);hold all;%figure(gcf);



