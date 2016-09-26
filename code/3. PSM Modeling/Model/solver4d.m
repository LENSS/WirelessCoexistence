% m=4;
% n=3;
clear all;
format long;

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

W_a=16; %contention window size
m_a=2; %the CW max.
l_a=2; %the retry limit.
s_a = 3; %ATIM packet size (slots)
w_a = 32; %ATIM window size (slots)

nn = 6; %number of total nodes, including myself
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


w0 = W; %contention window size (slots)
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

global atim_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global atim_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_atim_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_atim_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container

for i=1:l_a+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global A',num2str(i),';']);
    eval(['global B',num2str(i),';']);
end

global wait_atimwin_ends; %waiting for the end of the ATIM window
global wait_atimwin_ends2; %account for the extra states for other nodes' transmission during the waiting period (for atim window ends) of a node.
global wait_atimwin_ends_temp;
global wait_atimwin_ends_temp2;%this is for temporarily store the extra states for other nodes' transmission, see wait_atimwin_ends2.
global sum_wait_atimwin_ends; %for storing the sum of wait_atimwin_ends

if SYMBOL_VERIFICATION
    wait_atimwin_ends = sym(zeros((a-s+1)*nn,1)); %store the ATIM waiting chain
    wait_atimwin_ends2 = sym(zeros((a-s+1)*nn,1)); 
else
    wait_atimwin_ends = (zeros((a-s+1)*nn,1)); %store the ATIM waiting chain
    wait_atimwin_ends2 = (zeros((a-s+1)*nn,1)); 
end

%for DATA window
%=========================================


global data_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global data_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_data_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_data_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container

for i=1:l_d+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global C',num2str(i),';']);
    eval(['global D',num2str(i),';']);
end

global wait_datawin_ends;
global wait_datawin_ends_temp;
global wait_datawin_ends2;
global wait_datawin_ends_temp2;
global sum_wait_datawin_ends; %for storing the sum of wait_datawin_ends


if SYMBOL_VERIFICATION
    wait_datawin_ends = sym(zeros((w_d-s_d+1)*nn,1)); %store the DATA waiting chain
    wait_datawin_ends2 = sym(zeros((w_d-s_d+1)*nn,1)); 
else
    wait_datawin_ends = (zeros((w_d-s_d+1)*nn,1)); %store the DATA waiting chain
    wait_datawin_ends2 = (zeros((w_d-s_d+1)*nn,1)); 
end

    
for N=WIFI_START:5:WIFI_END

    %=================================================================
    %=================================================================
    %=================================================================
    %=================================================================
    %=================================================================
    %the following is the ATIM window part
    
    a=w_a;%20; %DATA window size (slots)
    s=s_a;%4; %DATA packet size (slots)
    W=W_a;
    l=l_a;
    m=m_a;

    
    
    for i=1:l+1
        if i<m+1
            if SYMBOL_VERIFICATION
                eval(['A',num2str(i),'=sym(zeros(a*nn,W*2^(i-1)));']);
                eval(['B',num2str(i),'=sym(zeros(s*nn+1,W*2^(i-1)));']);
            else
                eval(['A',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
                eval(['B',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
            end
        else
            if SYMBOL_VERIFICATION
                eval(['A',num2str(i),'=sym(zeros(a*nn,W*2^m));']);
                eval(['B',num2str(i),'=sym(zeros(s*nn+1,W*2^m));']);
            else
                eval(['A',num2str(i),'=zeros(a*nn,W*2^m);']);
                eval(['B',num2str(i),'=zeros(s*nn+1,W*2^m);']);
            end
        end
%         atim_backoff_states_container{i}=eval(['A',num2str(i)]);
%         atim_windowends_states_container{i}=eval(['B',num2str(i)]);
% %        disp(atim_backoff_states_container{1,i});
% %        disp(atim_windowends_states_container{1,i});
    end
    
    a=w_d;%20; %DATA window size (slots)
    s=s_d;%4; %DATA packet size (slots)
    W=W_d;
    l=l_d;
    m=m_d;
 
    
     for i=1:l+1
        if i<m+1
            eval(['C',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
            eval(['D',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
        else
            eval(['C',num2str(i),'=zeros(a*nn,W*2^m);']);
            eval(['D',num2str(i),'=zeros(s*nn+1,W*2^m);']);
        end
     end

    
    %     return;
    
%     if SYMBOL_VERIFICATION
% %         syms p_a p_b p_c p_d;
% 
%     else
%         tau=0.4;
%         tauc=1-tau;
%         
%         p_a=@(no) tauc^(no-1);
%         p_c=@(no) (no-1)*tau*tauc^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
%         p_b=@(no) 1-p_a(no)-p_c(no);
%         p_d=@(no) p_b(no)+p_c(no);
%     end
%     
% %     wait_atimwin_ends = 0;
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
%         end
%            
%         if i==1
%             for j=k:-1:1
%                 atim_backoff_states_container{1,i}(1,j) = 1/k;
%             end
%         else
%             Ct=atim_backoff_states_container{1,i-1}(:,1);
% 
% %                 disp(simplify(Ct));
% 
%             for j=k:-1:1
%                 for n=1:a*nn
%                     if Ct(n,1)~=0
%                         if mod(n,nn)==0
%                             no = 1;
%                         else
%                             no = nn - mod(n,nn) + 1;
%                         end
%                         if n+s*nn<=a*nn
%                             atim_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p_d(no)/k+atim_backoff_states_container{1,i}(n+s*nn,j);
%                         end
%                     end
%                 end
%             end
%         end
% % 
% %         disp('atim_backoff_states_container{1,i}');
% %         disp(atim_backoff_states_container{1,i});
% 
% %         return;
%         
%        for j=k:-1:1 %handle each column in the Markov Chain from right to left.
% 
%             for n=1:a*nn
%             %             disp(m);
%             %             disp(n);
%                 if atim_backoff_states_container{1,i}(n,j)~=0
%                     
%                     if mod(n,nn)==0
%                         no = 1;
%                     else
%                         no = nn - mod(n,nn) + 1;
%                     end
%                     
% %                     fprintf('no: %d\n', no);
%                     
%                     if j~=1
% %                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
%                         if n<=(a-s)*nn
% %                             fprintf('%s\n', char(atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no))));
%                             sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no));
%                             %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
%                         else
% %                             fprintf('(1+(a-n)*p_d(no) %d\n', 1+(a-ceil(n/nn)));
%                             sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*p_d(no));
%                             %here since if ATIM window ends, the freezing time may not be s-1 long.
%                         end
%                     else %j==1
% %                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
%                         if n<=(a-s)*nn
%                             sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1));
%                         else
%                             sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*1);
%                         end
%                     end
%                     
%                     if SYMBOL_VERIFICATION
% %                         disp('sum_atim_backoff_states:');
% %                         disp(simplify(sum_atim_backoff_states));
%                     else
% %                         disp('sum_atim_backoff_states:');
% %                         disp((sum_atim_backoff_states));
%                     end
% 
%                     if j>1 %except the first column
%             %                     disp('m>1');
%                         if n+1*nn<=a*nn
%                             atim_backoff_states_container{1,i}(n+1*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p_a(no)+atim_backoff_states_container{1,i}(n+1*nn,j-1);
%                         end          
%                         if n+s*nn<=a*nn
%                             atim_backoff_states_container{1,i}(n+s*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p_b(no)+atim_backoff_states_container{1,i}(n+s*nn,j-1);
%                         end
%                         if n+s*nn+1<=a*nn
%                             atim_backoff_states_container{1,i}(n+s*nn+1,j-1) = atim_backoff_states_container{1,i}(n,j)*p_c(no)+atim_backoff_states_container{1,i}(n+s*nn+1,j-1);
%                         end
%                     end
% 
%                     if n+1*nn>a*nn
% %                         fprintf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
%                         atim_windowends_states_container{1,i}(n+1*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p_a(no)+atim_windowends_states_container{1,i}(n+1*nn-a*nn,j);
%                     end
% 
%                     if n+s*nn>a*nn
% %                         fprintf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
%                         atim_windowends_states_container{1,i}(n+s*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p_b(no)+atim_windowends_states_container{1,i}(n+s*nn-a*nn,j);
%                     end
%                     
%                     if n+s*nn+1>a*nn
% %                         fprintf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
%                         atim_windowends_states_container{1,i}(n+s*nn-a*nn+1,j) = atim_backoff_states_container{1,i}(n,j)*p_c(no)+atim_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
%                     end
%                 end
%             end
%        end
%            
%             
% 
%         if SYMBOL_VERIFICATION
%             disp('atim_backoff_states_container{1,i} after...');
% 
%             disp(simplify(atim_backoff_states_container{1,i}));
% 
%             disp('atim_windowends_states_container{1,i} after...');
%             disp(simplify(atim_windowends_states_container{1,i}));
% 
%             disp('sum_atim_backoff_states:');
%             disp(simplify(sum_atim_backoff_states));
% 
%             wait_atimwin_ends_temp = sym(zeros((a-s+1)*nn,1)); %reset wait_atimwin_ends_temp
%             wait_atimwin_ends_temp2 = sym(zeros((a-s+1)*nn,1)); %reset wait_atimwin_ends_temp
%         else
%             disp('atim_backoff_states_container{1,i} after...');
%             disp((atim_backoff_states_container{1,i}));
% 
%             disp('atim_windowends_states_container{1,i} after...');
%             disp((atim_windowends_states_container{1,i}));
% 
%             disp('sum_atim_backoff_states:');
%             disp((sum_atim_backoff_states));
% 
%             wait_atimwin_ends_temp = zeros((a-s+1)*nn,1); %reset wait_atimwin_ends_temp
%             wait_atimwin_ends_temp2 = zeros((a-s+1)*nn,1); %reset wait_atimwin_ends_temp
%         end
%         
% %         return;
% 
%         for n=1:(a-s+1)*nn
%         %             disp(m);
%         %             disp(n);
%             if 1%atim_backoff_states_container{1,i}(n,1)~=0
%                 
%                 if mod(n,nn)==0
%                     no = 1;
%                 else
%                     no = nn - mod(n,nn) + 1;
%                 end
% 
%                 if n<=nn %one input
%                     wait_atimwin_ends_temp(n,1) = atim_backoff_states_container{1,i}(n,1)*p_a(no);
% %                     fprintf('%d:%d %g=%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),atim_backoff_states_container{1,i}(n,1),p_a(no));
%                 else
%                     if n<=nn*s %two inputs
%                         
%                         wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+ ...
%                             atim_backoff_states_container{1,i}(n,1)*p_a(no);
%                         
% %                         fprintf('%d:%d %g=%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
% %                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), ...
% %                             atim_backoff_states_container{1,i}(n,1),p_a(no));
%                     else
%                         if n==nn*s+1 % three inputs
%                             
%                         
%                             wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+...
%                                 wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)+ ...
%                                 atim_backoff_states_container{1,i}(n,1)*p_a(no);
%                             
%                             wait_atimwin_ends_temp2(n,1) = wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)*(s-1);
%                             
% %                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
% %                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), wait_atimwin_ends_temp(n-s*nn,1),p_b(no) ,...
% %                             atim_backoff_states_container{1,i}(n,1),p_a(no));
%                             
%                         else % four inputs
%                             
%                             if no+1==nn
%                                 tempno = no+1;
%                             else
%                                 tempno = mod(no+1,nn);
%                             end
%                             
%                             wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+ ...
%                                 wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)+ ...
%                                + wait_atimwin_ends_temp(n-s*nn-1,1)*p_c(tempno)+ ...
%                                atim_backoff_states_container{1,i}(n,1)*p_a(no);
%                            
%                             wait_atimwin_ends_temp2(n,1) = wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)*(s-1)+ ...
%                                + wait_atimwin_ends_temp(n-s*nn-1,1)*p_c(tempno)*(s-1);    
%                                
% %                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
% %                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), wait_atimwin_ends_temp(n-s*nn,1),p_b(no) ,...
% %                             wait_atimwin_ends_temp(n-s*nn-1,1),p_c(tempno), atim_backoff_states_container{1,i}(n,1),p_a(no));
%                         
%                         end
%                         
%                     end
%                 end
%             end
%         end
%         
%         
%          if SYMBOL_VERIFICATION
%         
% %             disp('wait_atimwin_ends_temp:');
% %             disp(simplify(wait_atimwin_ends_temp));
%          else
% %             disp('wait_atimwin_ends_temp:');
% %             disp((wait_atimwin_ends_temp));
%          end
% 
%         
%         wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
%         wait_atimwin_ends2 = wait_atimwin_ends2 + wait_atimwin_ends_temp2;
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
%    
%     disp('================ATIM summary================');
% 
%     
%     if SYMBOL_VERIFICATION
%         disp('wait_atimwin_ends:');
%         disp(simplify(wait_atimwin_ends));
%         
%         disp('wait_atimwin_ends2:');
%         disp(simplify(wait_atimwin_ends2));
%         
%         disp('sum_atim_backoff_states:');
%         sum_atim_backoff_states = simplify(sum_atim_backoff_states);
%         disp(sum_atim_backoff_states);   
%         
%         disp('sum_atim_windowends_states:');
%         sum_atim_windowends_states = simplify(sum_atim_windowends_states);
%         disp(sum_atim_windowends_states);
%         
%         sum_wait_atimwin_ends = simplify(sum(wait_atimwin_ends)+ sum(wait_atimwin_ends2));
%         disp('sum_wait_atimwin_ends:');
%         disp(sum_wait_atimwin_ends);
%     else
%         
%         disp('wait_atimwin_ends:');
%         disp((wait_atimwin_ends));
%         
%         disp('wait_atimwin_ends2:');
%         disp((wait_atimwin_ends2));
%             
%         disp('sum_atim_backoff_states:');
%         disp(sum_atim_backoff_states);   
%         
%         disp('sum_atim_windowends_states:');
%         disp(sum_atim_windowends_states);
%         
%         sum_wait_atimwin_ends = sum(wait_atimwin_ends) + sum(wait_atimwin_ends2);
%         disp('sum_wait_atimwin_ends:');
%         disp(sum_wait_atimwin_ends);
%         
%         all = sum_atim_backoff_states + sum_atim_windowends_states + sum_wait_atimwin_ends;
% %         fprintf('all: %g\n', all);
%         
% %         zw0 = 0;
% %         zw1 = 0;
% %         for i=1:l+1
% %             
% %             zw0 = zw0 + (sum(atim_backoff_states_container{1,i}(:,1)));
% %             zw1 = zw1 + (sum(sum(atim_backoff_states_container{1,i})));
% %             
% %             
% %         end
% %             tau = zw0/zw1;
% %             fprintf('tau: %g\n', tau);
%         
%     end
%     
%     
% %     return;
%     
%     %=================================================================
%     %=================================================================
%     %=================================================================
%     %=================================================================
%     %=================================================================
%     %the following is the DATA window part
%     %for test purpose, firstly use the same parameters as ATIM window
%     disp('******************************DATA start************************************');
% 
%     a=d;%20; %DATA window size (slots)
%     s=s_d;%4; %DATA packet size (slots)
%     W=W_d;
%     l=l_d;
%     m=m_d;
%     
%      for i=1:l+1
%         if i<m+1
%             if SYMBOL_VERIFICATION
%                 eval(['C',num2str(i),'=sym(zeros(a*nn,W*2^(i-1)));']);
%                 eval(['D',num2str(i),'=sym(zeros(s*nn+1,W*2^(i-1)));']);
%             else
%                 eval(['C',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
%                 eval(['D',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
%             end
%         else
%             if SYMBOL_VERIFICATION
%                 eval(['C',num2str(i),'=sym(zeros(a*nn,W*2^m));']);
%                 eval(['D',num2str(i),'=sym(zeros(s*nn+1,W*2^m));']);
%             else
%                 eval(['C',num2str(i),'=zeros(a*nn,W*2^m);']);
%                 eval(['D',num2str(i),'=zeros(s*nn+1,W*2^m);']);
%             end
%         end
%         data_backoff_states_container{i}=eval(['C',num2str(i)]);
%         data_windowends_states_container{i}=eval(['D',num2str(i)]);
% %        disp(data_backoff_states_container{1,i});
% %        disp(data_windowends_states_container{1,i});
%     end
%     
%     %     return;
%     
%     if SYMBOL_VERIFICATION
% %         syms p_a p_b p_c p_d;
% 
%     else
%         tau=0.2;
%         tauc=1-tau;
%         
%         p_a=@(no) tauc^(no-1);
%         p_c=@(no) (no-1)*tau*tauc^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
%         p_b=@(no) 1-p_a(no)-p_c(no);
%         p_d=@(no) p_b(no)+p_c(no);
%     end
%     
% %     wait_datawin_ends = 0;
%     sum_data_backoff_states = 0;
%     sum_data_windowends_states = 0;
%     sum_wait_datawin_ends = 0;
% 
%     for i=1:l+1
%         
%         if i<m+1
%             k = W*2^(i-1);%cw length
%          else %i>=m+1
%             k = W*2^m;%cw length
%         end
%            
%         if i==1
%             [temp_size temp]=size(wait_atimwin_ends);
%             for j=k:-1:1
% %                 fprintf('j: %d temp_size: %d\n', j, temp_size);
%                 for temp_j=1:nn
%                     data_backoff_states_container{1,i}(temp_j,j) = wait_atimwin_ends(temp_size-temp_j+1,1)/k; %%%%% TODO: 
%                 end
%             end
%         else
%             Ct=data_backoff_states_container{1,i-1}(:,1);
% 
% %                 disp(simplify(Ct));
% 
%             for j=k:-1:1
%                 for n=1:a*nn
%                     if Ct(n,1)~=0
%                         if mod(n,nn)==0
%                             no = 1;
%                         else
%                             no = nn - mod(n,nn) + 1;
%                         end
%                         if n+s*nn<=a*nn
%                             data_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p_d(no)/k+data_backoff_states_container{1,i}(n+s*nn,j);
%                         end
%                     end
%                 end
%             end
%         end
% % 
% %         disp('data_backoff_states_container{1,i}');
% %         disp(data_backoff_states_container{1,i});
% 
% %         return;
%         
%        for j=k:-1:1 %handle each column in the Markov Chain from right to left.
% 
%             for n=1:a*nn
%             %             disp(m);
%             %             disp(n);
%                 if data_backoff_states_container{1,i}(n,j)~=0
%                     
%                     if mod(n,nn)==0
%                         no = 1;
%                     else
%                         no = nn - mod(n,nn) + 1;
%                     end
%                     
% %                     fprintf('no: %d\n', no);
%                     
%                     if j~=1
% %                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
%                         if n<=(a-s)*nn
% %                             fprintf('%s\n', char(data_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no))));
%                             sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no));
%                             %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
%                         else
% %                             fprintf('(1+(a-n)*p_d(no) %d\n', 1+(a-ceil(n/nn)));
%                             sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*p_d(no));
%                             %here since if ATIM window ends, the freezing time may not be s-1 long.
%                         end
%                     else %j==1
% %                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
%                         if n<=(a-s)*nn
%                             sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(s-1));
%                         else
%                             sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*1);
%                         end
%                     end
%                     
%                     if SYMBOL_VERIFICATION
% %                         disp('sum_data_backoff_states:');
% %                         disp(simplify(sum_data_backoff_states));
%                     else
% %                         disp('sum_data_backoff_states:');
% %                         disp((sum_data_backoff_states));
%                     end
% 
%                     if j>1 %except the first column
%             %                     disp('m>1');
%                         if n+1*nn<=a*nn
%                             data_backoff_states_container{1,i}(n+1*nn,j-1) = data_backoff_states_container{1,i}(n,j)*p_a(no)+data_backoff_states_container{1,i}(n+1*nn,j-1);
%                         end          
%                         if n+s*nn<=a*nn
%                             data_backoff_states_container{1,i}(n+s*nn,j-1) = data_backoff_states_container{1,i}(n,j)*p_b(no)+data_backoff_states_container{1,i}(n+s*nn,j-1);
%                         end
%                         if n+s*nn+1<=a*nn
%                             data_backoff_states_container{1,i}(n+s*nn+1,j-1) = data_backoff_states_container{1,i}(n,j)*p_c(no)+data_backoff_states_container{1,i}(n+s*nn+1,j-1);
%                         end
%                     end
% 
%                     if n+1*nn>a*nn
% %                         fprintf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
%                         data_windowends_states_container{1,i}(n+1*nn-a*nn,j) = data_backoff_states_container{1,i}(n,j)*p_a(no)+data_windowends_states_container{1,i}(n+1*nn-a*nn,j);
%                     end
% 
%                     if n+s*nn>a*nn
% %                         fprintf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
%                         data_windowends_states_container{1,i}(n+s*nn-a*nn,j) = data_backoff_states_container{1,i}(n,j)*p_b(no)+data_windowends_states_container{1,i}(n+s*nn-a*nn,j);
%                     end
%                     
%                     if n+s*nn+1>a*nn
% %                         fprintf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
%                         data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j) = data_backoff_states_container{1,i}(n,j)*p_c(no)+data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
%                     end
%                 end
%             end
%        end
%            
%             
% 
%         if SYMBOL_VERIFICATION
%             disp('data_backoff_states_container{1,i} after...');
% 
%             disp(simplify(data_backoff_states_container{1,i}));
% 
%             disp('data_windowends_states_container{1,i} after...');
%             disp(simplify(data_windowends_states_container{1,i}));
% 
%             disp('sum_data_backoff_states:');
%             disp(simplify(sum_data_backoff_states));
% 
%             wait_datawin_ends_temp = sym(zeros((a-s+1)*nn,1)); %reset wait_datawin_ends_temp
%             wait_datawin_ends_temp2 = sym(zeros((a-s+1)*nn,1)); %reset wait_datawin_ends_temp
%         else
%             disp('data_backoff_states_container{1,i} after...');
%             disp((data_backoff_states_container{1,i}));
% 
%             disp('data_windowends_states_container{1,i} after...');
%             disp((data_windowends_states_container{1,i}));
% 
%             disp('sum_data_backoff_states:');
%             disp((sum_data_backoff_states));
% 
%             wait_datawin_ends_temp = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp
%             wait_datawin_ends_temp2 = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp
%         end
%         
% %         return;
% 
%         for n=1:(a-s+1)*nn
%         %             disp(m);
%         %             disp(n);
%             if 1%data_backoff_states_container{1,i}(n,1)~=0
%                 
%                 if mod(n,nn)==0
%                     no = 1;
%                 else
%                     no = nn - mod(n,nn) + 1;
%                 end
% 
%                 if n<=nn %one input
%                     wait_datawin_ends_temp(n,1) = data_backoff_states_container{1,i}(n,1)*p_a(no);
% %                     fprintf('%d:%d %g=%g*%g\n', n, no, wait_datawin_ends_temp(n,1),data_backoff_states_container{1,i}(n,1),p_a(no));
%                 else
%                     if n<=nn*s %two inputs
%                         
%                         wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+ ...
%                             data_backoff_states_container{1,i}(n,1)*p_a(no);
%                         
% %                         fprintf('%d:%d %g=%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
% %                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), ...
% %                             data_backoff_states_container{1,i}(n,1),p_a(no));
%                     else
%                         if n==nn*s+1 % three inputs
%                             
%                         
%                             wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+...
%                                 wait_datawin_ends_temp(n-s*nn,1)*p_b(no)+ ...
%                                 data_backoff_states_container{1,i}(n,1)*p_a(no);
%                             wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp(n-s*nn,1)*p_b(no)*(s-1);
% 
%                             
% %                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
% %                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), wait_datawin_ends_temp(n-s*nn,1),p_b(no) ,...
% %                             data_backoff_states_container{1,i}(n,1),p_a(no));
%                             
%                         else % four inputs
%                             
%                             if no+1==nn
%                                 tempno = no+1;
%                             else
%                                 tempno = mod(no+1,nn);
%                             end
%                             
%                             wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+ ...
%                                 wait_datawin_ends_temp(n-s*nn,1)*p_b(no)+ ...
%                                + wait_datawin_ends_temp(n-s*nn-1,1)*p_c(tempno)+ ...
%                                data_backoff_states_container{1,i}(n,1)*p_a(no);
%                             wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp(n-s*nn,1)*p_b(no)*(s-1)+ ...
%                                + wait_datawin_ends_temp(n-s*nn-1,1)*p_c(tempno)*(s-1);    
%                            
% %                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
% %                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), wait_datawin_ends_temp(n-s*nn,1),p_b(no) ,...
% %                             wait_datawin_ends_temp(n-s*nn-1,1),p_c(tempno), data_backoff_states_container{1,i}(n,1),p_a(no));
%                         
%                         end
%                         
%                     end
%                 end
%             end
%         end
%         
%         
%          if SYMBOL_VERIFICATION
%         
% %             disp('wait_datawin_ends_temp:');
% %             disp(simplify(wait_datawin_ends_temp));
%          else
% %             disp('wait_datawin_ends_temp:');
% %             disp((wait_datawin_ends_temp));
%          end
% 
%         
%         wait_datawin_ends = wait_datawin_ends + wait_datawin_ends_temp;
%         wait_datawin_ends2 = wait_datawin_ends2 + wait_datawin_ends_temp2;
%         
% %         disp('data_backoff_states_container:');
% %         disp(simplify(data_backoff_states_container{1,i}));
% %         disp('data_windowends_states_container:');
% %         disp(simplify(data_windowends_states_container{1,i}));
% %         
%         sum_data_windowends_states = sum_data_windowends_states + sum(sum(data_windowends_states_container{1,i}));
% %         disp('sum_data_windowends_states:');
% %         disp(sum_data_windowends_states);
% %         
% %         disp('wait_datawin_ends:');
% %         disp(simplify(wait_datawin_ends));
%     end
%     
%     disp('================DATA summary================');
% 
%     
%     if SYMBOL_VERIFICATION
%         disp('wait_datawin_ends:');
%         disp(simplify(wait_datawin_ends));
%         
%         disp('wait_datawin_ends2:');
%         disp(simplify(wait_datawin_ends2));
%         
%         disp('sum_data_backoff_states:');
%         sum_data_backoff_states = simplify(sum_data_backoff_states);
%         disp(sum_data_backoff_states);   
%         
%         disp('sum_data_windowends_states:');
%         sum_data_windowends_states = simplify(sum_data_windowends_states);
%         disp(sum_data_windowends_states);
%         
%         sum_wait_datawin_ends = simplify(sum(wait_datawin_ends)+sum(wait_datawin_ends2));
%         disp('sum_wait_datawin_ends:');
%         disp(sum_wait_datawin_ends);
%     else
%         
%         disp('wait_datawin_ends:');
%         disp((wait_datawin_ends));
%         
%         disp('wait_datawin_ends2:');
%         disp((wait_datawin_ends2));
%             
%         disp('sum_data_backoff_states:');
%         disp(sum_data_backoff_states);   
%         
%         disp('sum_data_windowends_states:');
%         disp(sum_data_windowends_states);
%         
%         sum_wait_datawin_ends = sum(wait_datawin_ends) + sum(wait_datawin_ends2);
%         disp('sum_wait_datawin_ends:');
%         disp(sum_wait_datawin_ends);
%         
%         all = all + sum_data_backoff_states + sum_data_windowends_states + sum_wait_datawin_ends;
%         fprintf('all: %g\n', all);
%         
% %         zw0 = 0;
% %         zw1 = 0;
% %         for i=1:l+1
% %             
% %             zw0 = zw0 + (sum(data_backoff_states_container{1,i}(:,1)));
% %             zw1 = zw1 + (sum(sum(data_backoff_states_container{1,i})));
% %             
% %         end
% %             tau = zw0/zw1;
% %             fprintf('tau: %g\n', tau);
% 
%     end   
%     
%         
%     return;
% 

    
    xguess=[0.2 0.2 0.1]';
    options = optimset('MaxFunEvals',2000, 'MaxIter', 2000);
    xvect = fsolve('solver4d_math', xguess, options);
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



