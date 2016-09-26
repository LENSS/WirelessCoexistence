
function f=solver4d_math(xvect)
global N;
% global Ls;
global Lo;
global W;
global m;
global l;
global s;
global a;
global nn;

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

global L_s;
global L_c;
global L_sc;

global atim_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global atim_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_atim_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_atim_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container
global wait_atimwin_ends; %waiting for the end of the ATIM window
global wait_atimwin_ends2; %waiting for the end of the ATIM window
global wait_atimwin_ends_temp;
global wait_atimwin_ends_temp2;
global sum_wait_atimwin_ends; %for storing the sum of wait_atimwin_ends

global data_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global data_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_data_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_data_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container
global wait_datawin_ends;
global wait_datawin_ends_temp;
global wait_datawin_ends2;
global wait_datawin_ends_temp2;
global sum_wait_datawin_ends; %for storing the sum of wait_datawin_ends


for i=1:l_a+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global A',num2str(i),';']);
    eval(['global B',num2str(i),';']);
%     eval(['A',num2str(i),'=zeros(1,W*2^(i-1))']);
%     B{i}=eval(['A',num2str(i)]);
end
for i=1:l_d+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global C',num2str(i),';']);
    eval(['global D',num2str(i),';']);
end

    tau = xvect(1);
    tau2 = xvect(2);
    init = xvect(3);
%     tau = xvect(3);

    a=w_a;%20; %DATA window size (slots)
    s=s_a;%4; %DATA packet size (slots)
    W=W_a;
    l=l_a;
    m=m_a;

    for i=1:l+1
        atim_backoff_states_container{i}=eval(['A',num2str(i)]);
        atim_windowends_states_container{i}=eval(['B',num2str(i)]);
    end
    
%     p=0.2;
%     wait_atimwin_ends = 0;
    wait_atimwin_ends = (zeros((a-s+1)*nn,1)); %store the ATIM waiting chain
    wait_atimwin_ends2 = (zeros((a-s+1)*nn,1)); 
    sum_atim_backoff_states = 0;
    sum_atim_windowends_states = 0;
    sum_wait_atimwin_ends = 0;
    
    
    tauc=1-tau;

    p_a=@(no) tauc^(no-1);
    p_c=@(no) (no-1)*tau*tauc^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
    p_b=@(no) 1-p_a(no)-p_c(no);
    p_d=@(no) p_b(no)+p_c(no);
    
%     wait_atimwin_ends = 0;
    sum_atim_backoff_states = 0;
    sum_atim_windowends_states = 0;
    sum_wait_atimwin_ends = 0;

    for i=1:l+1
        
        if i<m+1
            k = W*2^(i-1);%cw length
         else %i>=m+1
            k = W*2^m;%cw length
        end
           
        if i==1
            for j=k:-1:1
                atim_backoff_states_container{1,i}(1,j) = 1/k;
            end
        else
            Ct=atim_backoff_states_container{1,i-1}(:,1);

%                 disp(simplify(Ct));

            for j=k:-1:1
                for n=1:a*nn
                    if Ct(n,1)~=0
                        if mod(n,nn)==0
                            no = 1;
                        else
                            no = nn - mod(n,nn) + 1;
                        end
                        if n+s*nn<=a*nn
                            atim_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p_d(no)/k+atim_backoff_states_container{1,i}(n+s*nn,j);
                        end
                    end
                end
            end
        end
% 
%         disp('atim_backoff_states_container{1,i}');
%         disp(atim_backoff_states_container{1,i});

%         return;
        
       for j=k:-1:1 %handle each column in the Markov Chain from right to left.

            for n=1:a*nn
            %             disp(m);
            %             disp(n);
                if atim_backoff_states_container{1,i}(n,j)~=0
                    
                    if mod(n,nn)==0
                        no = 1;
                    else
                        no = nn - mod(n,nn) + 1;
                    end
                    
%                     fprintf('no: %d\n', no);
                    
                    if j~=1
%                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
                        if n<=(a-s)*nn
%                             fprintf('%s\n', char(atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no))));
                            sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no));
                            %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
                        else
%                             fprintf('(1+(a-n)*p_d(no) %d\n', 1+(a-ceil(n/nn)));
                            sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*p_d(no));
                            %here since if ATIM window ends, the freezing time may not be s-1 long.
                        end
                    else %j==1
%                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
                        if n<=(a-s)*nn
                            sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1));
                        else
                            sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*1);
                        end
                    end
                    

%                         disp('sum_atim_backoff_states:');
%                         disp((sum_atim_backoff_states));

                    if j>1 %except the first column
            %                     disp('m>1');
                        if n+1*nn<=a*nn
                            atim_backoff_states_container{1,i}(n+1*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p_a(no)+atim_backoff_states_container{1,i}(n+1*nn,j-1);
                        end          
                        if n+s*nn<=a*nn
                            atim_backoff_states_container{1,i}(n+s*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p_b(no)+atim_backoff_states_container{1,i}(n+s*nn,j-1);
                        end
                        if n+s*nn+1<=a*nn
                            atim_backoff_states_container{1,i}(n+s*nn+1,j-1) = atim_backoff_states_container{1,i}(n,j)*p_c(no)+atim_backoff_states_container{1,i}(n+s*nn+1,j-1);
                        end
                    end

                    if n+1*nn>a*nn
%                         fprintf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
                        atim_windowends_states_container{1,i}(n+1*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p_a(no)+atim_windowends_states_container{1,i}(n+1*nn-a*nn,j);
                    end

                    if n+s*nn>a*nn
%                         fprintf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
                        atim_windowends_states_container{1,i}(n+s*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p_b(no)+atim_windowends_states_container{1,i}(n+s*nn-a*nn,j);
                    end
                    
                    if n+s*nn+1>a*nn
%                         fprintf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
                        atim_windowends_states_container{1,i}(n+s*nn-a*nn+1,j) = atim_backoff_states_container{1,i}(n,j)*p_c(no)+atim_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
                    end
                end
            end
       end
           
            
% 
%             disp('atim_backoff_states_container{1,i} after...');
%             disp((atim_backoff_states_container{1,i}));
% 
%             disp('atim_windowends_states_container{1,i} after...');
%             disp((atim_windowends_states_container{1,i}));
% 
%             disp('sum_atim_backoff_states:');
%             disp((sum_atim_backoff_states));

            wait_atimwin_ends_temp = zeros((a-s+1)*nn,1); %reset wait_atimwin_ends_temp
            wait_atimwin_ends_temp2 = zeros((a-s+1)*nn,1); %reset wait_atimwin_ends_temp
        
%         return;

        for n=1:(a-s+1)*nn
        %             disp(m);
        %             disp(n);
            if 1%atim_backoff_states_container{1,i}(n,1)~=0
                
                if mod(n,nn)==0
                    no = 1;
                else
                    no = nn - mod(n,nn) + 1;
                end

                if n<=nn %one input
                    wait_atimwin_ends_temp(n,1) = atim_backoff_states_container{1,i}(n,1)*p_a(no);
%                     fprintf('%d:%d %g=%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),atim_backoff_states_container{1,i}(n,1),p_a(no));
                else
                    if n<=nn*s %two inputs
                        
                        wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+ ...
                            atim_backoff_states_container{1,i}(n,1)*p_a(no);
                        
%                         fprintf('%d:%d %g=%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
%                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), ...
%                             atim_backoff_states_container{1,i}(n,1),p_a(no));
                    else
                        if n==nn*s+1 % three inputs
                            
                        
                            wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+...
                                wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)+ ...
                                atim_backoff_states_container{1,i}(n,1)*p_a(no);
                            
                            wait_atimwin_ends_temp2(n,1) = wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)*(s-1);
                            
%                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
%                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), wait_atimwin_ends_temp(n-s*nn,1),p_b(no) ,...
%                             atim_backoff_states_container{1,i}(n,1),p_a(no));
                            
                        else % four inputs
                            
                            if no+1==nn
                                tempno = no+1;
                            else
                                tempno = mod(no+1,nn);
                            end
                            
                            wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1*nn,1)*p_a(no)+ ...
                                wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)+ ...
                               + wait_atimwin_ends_temp(n-s*nn-1,1)*p_c(tempno)+ ...
                               atim_backoff_states_container{1,i}(n,1)*p_a(no);
                           
                            wait_atimwin_ends_temp2(n,1) = wait_atimwin_ends_temp(n-s*nn,1)*p_b(no)*(s-1)+ ...
                               + wait_atimwin_ends_temp(n-s*nn-1,1)*p_c(tempno)*(s-1);    
                               
%                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g+%g*%g\n', n, no, wait_atimwin_ends_temp(n,1),...
%                             wait_atimwin_ends_temp(n-1*nn,1),p_a(no), wait_atimwin_ends_temp(n-s*nn,1),p_b(no) ,...
%                             wait_atimwin_ends_temp(n-s*nn-1,1),p_c(tempno), atim_backoff_states_container{1,i}(n,1),p_a(no));
                        
                        end
                        
                    end
                end
            end
        end
        
        
%             disp('wait_atimwin_ends_temp:');
%             disp((wait_atimwin_ends_temp));

        
        wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
        wait_atimwin_ends2 = wait_atimwin_ends2 + wait_atimwin_ends_temp2;
       
%         disp('atim_backoff_states_container:');
%         disp(simplify(atim_backoff_states_container{1,i}));
%         disp('atim_windowends_states_container:');
%         disp(simplify(atim_windowends_states_container{1,i}));
%         
%         
        sum_atim_windowends_states = sum_atim_windowends_states + sum(sum(atim_windowends_states_container{1,i}));
%         disp('sum_atim_windowends_states:');
%         disp(sum_atim_windowends_states);
%         
%         disp('wait_atimwin_ends:');
%         disp(simplify(wait_atimwin_ends));
    end
    
%    
%     disp('================ATIM summary================');
% 
%     
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
        
        sum_wait_atimwin_ends = sum(wait_atimwin_ends) + sum(wait_atimwin_ends2);
%         disp('sum_wait_atimwin_ends:');
%         disp(sum_wait_atimwin_ends);
        

        
        
%         all = sum_atim_backoff_states + sum_atim_windowends_states + sum_wait_atimwin_ends;
%         fprintf('all: %g\n', all);
        
        zw0 = 0;
        zw1 = 0;
        for i=1:l+1
            
            zw0 = zw0 + (sum(atim_backoff_states_container{1,i}(:,1)));
            zw1 = zw1 + (sum(sum(atim_backoff_states_container{1,i})));
            
            
        end
%         tau = zw0/zw1;
%         fprintf('tau: %g, init: %g, p: %g\n', tau, init, p);


    a=w_d;%20; %DATA window size (slots)
    s=s_d;%4; %DATA packet size (slots)
    W=W_d;
    l=l_d;
    m=m_d;
    
     for i=1:l+1
%         if i<m+1
%             eval(['C',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
%             eval(['D',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
%         else
%             eval(['C',num2str(i),'=zeros(a*nn,W*2^m);']);
%             eval(['D',num2str(i),'=zeros(s*nn+1,W*2^m);']);
%         end
        data_backoff_states_container{i}=eval(['C',num2str(i)]);
        data_windowends_states_container{i}=eval(['D',num2str(i)]);
%        disp(data_backoff_states_container{1,i});
%        disp(data_windowends_states_container{1,i});
    end
    
    %     return;
    
        tauc2=1-tau2;
        
        p_a=@(no) tauc2^(no-1);
        p_c=@(no) (no-1)*tau2*tauc2^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
        p_b=@(no) 1-p_a(no)-p_c(no);
        p_d=@(no) p_b(no)+p_c(no);
    
%     wait_datawin_ends = 0;
    wait_datawin_ends = (zeros((w_d-s_d+1)*nn,1)); %store the DATA waiting chain
    wait_datawin_ends2 = (zeros((w_d-s_d+1)*nn,1)); 
    sum_data_backoff_states = 0;
    sum_data_windowends_states = 0;
    sum_wait_datawin_ends = 0;

    for i=1:l+1
        
        if i<m+1
            k = W*2^(i-1);%cw length
         else %i>=m+1
            k = W*2^m;%cw length
        end
           
        if i==1
            [temp_size temp]=size(wait_atimwin_ends);
            for j=k:-1:1
%                 fprintf('j: %d temp_size: %d\n', j, temp_size);
                for temp_j=1:nn
                    data_backoff_states_container{1,i}(temp_j,j) = wait_atimwin_ends(temp_size-temp_j+1,1)/k; 
                end
            end
        else
            Ct=data_backoff_states_container{1,i-1}(:,1);

%                 disp(simplify(Ct));

            for j=k:-1:1
                for n=1:a*nn
                    if Ct(n,1)~=0
                        if mod(n,nn)==0
                            no = 1;
                        else
                            no = nn - mod(n,nn) + 1;
                        end
                        if n+s*nn<=a*nn
                            data_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p_d(no)/k+data_backoff_states_container{1,i}(n+s*nn,j);
                        end
                    end
                end
            end
        end
% 
%         disp('data_backoff_states_container{1,i}');
%         disp(data_backoff_states_container{1,i});

%         return;
        
       for j=k:-1:1 %handle each column in the Markov Chain from right to left.

            for n=1:a*nn
            %             disp(m);
            %             disp(n);
                if data_backoff_states_container{1,i}(n,j)~=0
                    
                    if mod(n,nn)==0
                        no = 1;
                    else
                        no = nn - mod(n,nn) + 1;
                    end
                    
%                     fprintf('no: %d\n', no);
                    
                    if j~=1
%                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
                        if n<=(a-s)*nn
%                             fprintf('%s\n', char(data_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no))));
                            sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no));
                            %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
                        else
%                             fprintf('(1+(a-n)*p_d(no) %d\n', 1+(a-ceil(n/nn)));
                            sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*p_d(no));
                            %here since if ATIM window ends, the freezing time may not be s-1 long.
                        end
                    else %j==1
%                             fprintf('i: %d, j: %d, n: %d, (a-s)*nn: %d\n', i,j, n, (a-s)*nn);
                        if n<=(a-s)*nn
                            sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(s-1));
                        else
                            sum_data_backoff_states = sum_data_backoff_states + data_backoff_states_container{1,i}(n,j)*(1+(a-ceil(n/nn))*1);
                        end
                    end
                    
%                         disp('sum_data_backoff_states:');
%                         disp((sum_data_backoff_states));

                    if j>1 %except the first column
            %                     disp('m>1');
                        if n+1*nn<=a*nn
                            data_backoff_states_container{1,i}(n+1*nn,j-1) = data_backoff_states_container{1,i}(n,j)*p_a(no)+data_backoff_states_container{1,i}(n+1*nn,j-1);
                        end          
                        if n+s*nn<=a*nn
                            data_backoff_states_container{1,i}(n+s*nn,j-1) = data_backoff_states_container{1,i}(n,j)*p_b(no)+data_backoff_states_container{1,i}(n+s*nn,j-1);
                        end
                        if n+s*nn+1<=a*nn
                            data_backoff_states_container{1,i}(n+s*nn+1,j-1) = data_backoff_states_container{1,i}(n,j)*p_c(no)+data_backoff_states_container{1,i}(n+s*nn+1,j-1);
                        end
                    end

                    if n+1*nn>a*nn
%                         fprintf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
                        data_windowends_states_container{1,i}(n+1*nn-a*nn,j) = data_backoff_states_container{1,i}(n,j)*p_a(no)+data_windowends_states_container{1,i}(n+1*nn-a*nn,j);
                    end

                    if n+s*nn>a*nn
%                         fprintf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
                        data_windowends_states_container{1,i}(n+s*nn-a*nn,j) = data_backoff_states_container{1,i}(n,j)*p_b(no)+data_windowends_states_container{1,i}(n+s*nn-a*nn,j);
                    end
                    
                    if n+s*nn+1>a*nn
%                         fprintf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
                        data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j) = data_backoff_states_container{1,i}(n,j)*p_c(no)+data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
                    end
                end
            end
       end
           
            
% 
%             disp('data_backoff_states_container{1,i} after...');
%             disp((data_backoff_states_container{1,i}));
% 
%             disp('data_windowends_states_container{1,i} after...');
%             disp((data_windowends_states_container{1,i}));
% 
%             disp('sum_data_backoff_states:');
%             disp((sum_data_backoff_states));

            wait_datawin_ends_temp = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp
            wait_datawin_ends_temp2 = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp
        
%         return;

        for n=1:(a-s+1)*nn
        %             disp(m);
        %             disp(n);
            if 1%data_backoff_states_container{1,i}(n,1)~=0
                
                if mod(n,nn)==0
                    no = 1;
                else
                    no = nn - mod(n,nn) + 1;
                end

                if n<=nn %one input
                    wait_datawin_ends_temp(n,1) = data_backoff_states_container{1,i}(n,1)*p_a(no);
%                     fprintf('%d:%d %g=%g*%g\n', n, no, wait_datawin_ends_temp(n,1),data_backoff_states_container{1,i}(n,1),p_a(no));
                else
                    if n<=nn*s %two inputs
                        
                        wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+ ...
                            data_backoff_states_container{1,i}(n,1)*p_a(no);
                        
%                         fprintf('%d:%d %g=%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
%                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), ...
%                             data_backoff_states_container{1,i}(n,1),p_a(no));
                    else
                        if n==nn*s+1 % three inputs
                            
                        
                            wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+...
                                wait_datawin_ends_temp(n-s*nn,1)*p_b(no)+ ...
                                data_backoff_states_container{1,i}(n,1)*p_a(no);
                            wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp(n-s*nn,1)*p_b(no)*(s-1);

                            
%                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
%                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), wait_datawin_ends_temp(n-s*nn,1),p_b(no) ,...
%                             data_backoff_states_container{1,i}(n,1),p_a(no));
                            
                        else % four inputs
                            
                            if no+1==nn
                                tempno = no+1;
                            else
                                tempno = mod(no+1,nn);
                            end
                            
                            wait_datawin_ends_temp(n,1) = wait_datawin_ends_temp(n-1*nn,1)*p_a(no)+ ...
                                wait_datawin_ends_temp(n-s*nn,1)*p_b(no)+ ...
                               + wait_datawin_ends_temp(n-s*nn-1,1)*p_c(tempno)+ ...
                               data_backoff_states_container{1,i}(n,1)*p_a(no);
                            wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp(n-s*nn,1)*p_b(no)*(s-1)+ ...
                               + wait_datawin_ends_temp(n-s*nn-1,1)*p_c(tempno)*(s-1);    
                           
%                         fprintf('%d:%d %g=%g*%g+%g*%g+%g*%g+%g*%g\n', n, no, wait_datawin_ends_temp(n,1),...
%                             wait_datawin_ends_temp(n-1*nn,1),p_a(no), wait_datawin_ends_temp(n-s*nn,1),p_b(no) ,...
%                             wait_datawin_ends_temp(n-s*nn-1,1),p_c(tempno), data_backoff_states_container{1,i}(n,1),p_a(no));
                        
                        end
                        
                    end
                end
            end
        end
        
        
%             disp('wait_datawin_ends_temp:');
%             disp((wait_datawin_ends_temp));

        
        wait_datawin_ends = wait_datawin_ends + wait_datawin_ends_temp;
        wait_datawin_ends2 = wait_datawin_ends2 + wait_datawin_ends_temp2;
        
%         disp('data_backoff_states_container:');
%         disp(simplify(data_backoff_states_container{1,i}));
%         disp('data_windowends_states_container:');
%         disp(simplify(data_windowends_states_container{1,i}));
%         
        sum_data_windowends_states = sum_data_windowends_states + sum(sum(data_windowends_states_container{1,i}));
%         disp('sum_data_windowends_states:');
%         disp(sum_data_windowends_states);
%         
%         disp('wait_datawin_ends:');
%         disp(simplify(wait_datawin_ends));
    end
    
%     disp('================DATA summary================');
% 
%     
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
        
        sum_wait_datawin_ends = sum(wait_datawin_ends) + sum(wait_datawin_ends2);
%         disp('sum_wait_datawin_ends:');
%         disp(sum_wait_datawin_ends);
        
        all = sum_atim_backoff_states + sum_atim_windowends_states*a + sum_wait_atimwin_ends + sum_data_backoff_states + sum_data_windowends_states + sum_wait_datawin_ends;
        fprintf('all: %g\n', all);


        zw3 = 0;
        zw4 = 0;
        for i=1:l+1
            
            zw3 = zw3 + (sum(data_backoff_states_container{1,i}(:,1)));
            zw4 = zw4 + (sum(sum(data_backoff_states_container{1,i})));
            
            
        end




    
%     f(1) =  1 - (1-tau)^(N-1) - p;
    f(1) =  zw0/zw1 - tau;
    f(2) =  zw3/zw4 - tau2;

    f(3) = 1/all - init;
    












