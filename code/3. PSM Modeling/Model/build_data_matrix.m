function [ output_args ] = build_data_matrix( zw )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global nn;
    global W_d;
    global m_d;
    global l_d;
    global s_d;
    global w_d;
    global w_a;
    global data_p_save_container;
    global data_pe_container;
    global data_pf_container;
    global data_ps_container;
    global throughput_container_temp;
    global data_energy_container_temp;
    global data_sumstate_container;
    
    a=w_d;
    s=s_d;
    W=W_d;
    l=l_d;
    m=m_d;

%     global data_backoff_states_container;% data_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
%     global data_windowends_states_container;% data_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
%     global sum_data_backoff_states;% sum_data_backoff_states is for storing the sum of data_backoff_states_container
%     global sum_data_windowends_states;% sum_data_backoff_states is for storing the sum of data_windowends_states_container
%     global wait_datawin_ends; %waiting for the end of the ATIM window
%     global wait_datawin_ends2; %waiting for the end of the ATIM window
%     global wait_datawin_ends_temp;
%     global wait_datawin_ends_temp2;
%     global sum_wait_datawin_ends; %for storing the sum of wait_datawin_ends

    global data_backoff_states_container;% data_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
    global data_windowends_states_container;% data_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
%     global sum_data_backoff_states;% sum_data_backoff_states is for storing the sum of data_backoff_states_container
%     global sum_data_windowends_states;% sum_data_backoff_states is for storing the sum of data_windowends_states_container
%     global wait_datawin_ends;
%     global wait_datawin_ends_temp;
%     global wait_datawin_ends2;
%     global wait_datawin_ends_temp2;
%     global sum_wait_datawin_ends; %for storing the sum of wait_datawin_ends

    sum_data_backoff_states = 0;% sum_data_backoff_states is for storing the sum of data_backoff_states_container
	sum_data_windowends_states = 0;% sum_data_backoff_states is for storing the sum of data_windowends_states_container
	sum_wait_datawin_ends = 0; %for storing the sum of wait_datawin_ends



%     global P;
    global PRECISION_INIT;
    global SYMBOL_VERIFICATION;
%     global DEBUG;
    
    global tau;
    global tauc;
%     global p1 p2 p3 p4;

for i=1:l_d+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global C',num2str(i),';']);
    eval(['global D',num2str(i),';']);
end
    
    p2e = zeros(nn,nn); 
    p2f = zeros(nn,nn); 
    p2s = zeros(nn,nn);
    init = 1;


%     pex = zeros(nn,nn); 
%     pnx = zeros(nn,nn); 
%     pwx = zeros(nn,nn); 


%     eval(['init',num2str(zw),'=0.1;']);
%     eval(['taua',num2str(zw),'=0.2;']);

    PRECISION = PRECISION_INIT*10^(-floor((nn-zw)/5));

    for i=1:l+1
        data_backoff_states_container{i}=eval(['C',num2str(i)]);
        data_windowends_states_container{i}=eval(['D',num2str(i)]);
    end    



    if SYMBOL_VERIFICATION
        wait_datawin_ends = sym(zeros((a-s+1)*nn,1)); %store the ATIM waiting chain
        wait_datawin_ends2 = sym(zeros((a-s+1)*nn,1)); 
    else
        wait_datawin_ends = (zeros((a-s+1)*nn,1)); %store the ATIM waiting chain
        wait_datawin_ends2 = (zeros((a-s+1)*nn,1)); 
    end


%     a=w_d;%20; %DATA window size (slots)
%     s=s_d;%4; %DATA packet size (slots)
%     W=W_d;
%     l=l_d;
%     m=m_d;
%  
%     
%      for i=1:l+1
%         if i<m+1
%             eval(['C',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
%             eval(['D',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
%         else
%             eval(['C',num2str(i),'=zeros(a*nn,W*2^m);']);
%             eval(['D',num2str(i),'=zeros(s*nn+1,W*2^m);']);
%         end
%      end


    %     return;
    NEWMC=3;
    %0. the original method to compute tau, slow and is not reusable (1024*10000*20 needs 1800s to get one matrix)
    %1. the most accurate way to compute tau, slow but can be reused (overall speed is the fastest, but to build the matrix is slow)
    %2. the modified version of the original method, faster but still needs some time (100x faster, but still needs 18s)
    %3. modified version of the most accurate way to compute tau

    if SYMBOL_VERIFICATION
%         syms p_a p_b p_c p_d;

    else

    end

%     wait_datawin_ends = 0;


    tt=tic;
    for i=1:l+1%This for loop is not used, after the first iteration, it will break

        if i<m+1
            k = W*2^(i-1);%cw length
         else %i>=m+1
            k = W*2^m;%cw length
        end

        if i==1
%             for j=k:-1:1
%                 init = eval(['init',num2str(zw),';']);
                data_backoff_states_container{1,i}(zw,1:k) = init/k;
%                     data_backoff_states_container{1,i}(zw,k) = init/k;
%             end
        else
            Ct=data_backoff_states_container{1,i-1}(:,1);
%                 zw_disp((data_backoff_states_container{1,i}));
            [row]=find(Ct>PRECISION);

%                 zw_disp(simplify(Ct));

%                 for j=k:-1:1
%                 for n=1:a*nn
%                     if Ct(n,1)~=0
                      for n=row.'
                        if mod(n,nn)==0
                            no = 1;
                        else
                            no = nn - mod(n,nn) + 1;
                        end
                        if n+s*nn<=a*nn
                            data_backoff_states_container{1,i}(n+s*nn,1:k) = Ct(n,1)*p_d(no)/k;% + ...
%                                     data_backoff_states_container{1,i}(n+s*nn,1:k);
                        end
                    end
%                 end
%                 end
        end
        wait_datawin_ends_temp = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp

%         zw_disp('data_backoff_states_container{1,i}');
%         zw_disp(data_backoff_states_container{1,i});

%         return;
%             tic;



        t2=tic;
        p_save = zeros(a*nn,4);


       for n=1:a*nn %handle each row in the Markov Chain from top to down.
%             for n=1:a*nn %handle every element in a column
%             %             zw_disp(m);
%             %             zw_disp(n);
%                 if data_backoff_states_container{1,i}(n,j)~=0
            if mod(n,nn)==0
                no = 1;
            else
                no = nn - mod(n,nn) + 1;
            end
            truetau = 0;
            denominator = 0;
            for j=1:l+1
                truetau=truetau+data_backoff_states_container{1,j}(n,1);
                denominator=denominator+sum(data_backoff_states_container{1,j}(n,:));
            end
            if truetau==0 || denominator==0
                continue;
            end
            tau = truetau/denominator;
%                 tau = 0.1;
            tauc=1-tau;

%                 p_a=@(no) tauc^(no-1);
%                 p_c=@(no) (no-1)*tau*tauc^(no-2);%nchoosek(no,1)*tau*tauc^(no-1);
%                 p_b=@(no) 1-p_a(no)-p_c(no);
%                 p_d=@(no) p_b(no)+p_c(no);

%                 for j=1:nn
               p_a(no)=p_af(no);
               p_b(no)=p_bf(no);
               p_c(no)=p_cf(no);
               p_d(no)=p_df(no);
%                 end
               p_save(n,1) = p_a(no);
               p_save(n,2) = p_b(no);
               p_save(n,3) = p_c(no);
               p_save(n,4) = p_d(no);

% %                 for j=1:nn
%                    p_a(j)=p_af(j);
%                    p_b(j)=p_bf(j);
%                    p_c(j)=p_cf(j);
%                    p_d(j)=p_df(j);
% %                 end

            t3=tic;
%                 zw_printf('[ %d ] tau: %f, ', n, tau);
            for j=1:l+1 %This loop is what is used in real.

                if j<m+1
                    k = W*2^(j-1);%cw length
                else %i>=m+1
                    k = W*2^m;%cw length
                end



                if n<=(a-s)*nn
%                             zw_printf('%s\n', char(data_backoff_states_container{1,i}(n,j)*(1+(s-1)*p_d(no))));
                    sum_data_backoff_states = sum_data_backoff_states + ...
                        sum(data_backoff_states_container{1,j}(n,2:k))*(1+(s-1)*p_d(no));
                    %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
                    sum_data_backoff_states = sum_data_backoff_states + ...
                        data_backoff_states_container{1,j}(n,1)*(1+(s-1));
                else
%                             zw_printf('(1+(a-n)*p_d(no) %d\n', 1+(a-ceil(n/nn)));
                    sum_data_backoff_states = sum_data_backoff_states + ...
                        sum(data_backoff_states_container{1,j}(n,2:k))*(1+(a-ceil(n/nn))*p_d(no));
                    %here since if ATIM window ends, the freezing time may not be s-1 long.
                    sum_data_backoff_states = sum_data_backoff_states + ...
                        data_backoff_states_container{1,j}(n,1)*(1+(a-ceil(n/nn))*1);
                end

                if SYMBOL_VERIFICATION
%                         zw_disp('sum_data_backoff_states:');
%                         zw_disp(simplify(sum_data_backoff_states));
                else
%                         zw_disp('sum_data_backoff_states:');
%                         zw_disp((sum_data_backoff_states));
                end

%                     if j>1 %except the first column
        %                     zw_disp('m>1');
                if n+1*nn<=a*nn
                    data_backoff_states_container{1,j}(n+1*nn,1:k-1) = data_backoff_states_container{1,j}(n,2:k)*p_a(no)+...
                        data_backoff_states_container{1,j}(n+1*nn,1:k-1);
                end          
                if n+s*nn<=a*nn
                    data_backoff_states_container{1,j}(n+s*nn,1:k-1) = data_backoff_states_container{1,j}(n,2:k)*p_b(no)+...
                        data_backoff_states_container{1,j}(n+s*nn,1:k-1);
                end
                if n+s*nn+1<=a*nn
                    data_backoff_states_container{1,j}(n+s*nn+1,1:k-1) = data_backoff_states_container{1,j}(n,2:k)*p_c(no)+...
                        data_backoff_states_container{1,j}(n+s*nn+1,1:k-1);
                end

                if n+1*nn>a*nn
%                         zw_printf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
%                                     temp=data_windowends_states_container{1,j}(n+1*nn-a*nn,j);
                    data_windowends_states_container{1,j}(n+1*nn-a*nn,2:k) = data_backoff_states_container{1,j}(n,2:k)*p_a(no)+...
                        data_windowends_states_container{1,j}(n+1*nn-a*nn,2:k);
%                         zw_disp(data_windowends_states_container{1,j}(n+1*nn-a*nn,2:k));
%                                     zw_printf('11. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+1*nn-a*nn,j, ...
%                                         data_windowends_states_container{1,i}(n+1*nn-a*nn,j),data_backoff_states_container{1,i}(n,j),p_a(no), temp);
                end

                if n+s*nn>a*nn
%                         zw_printf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
%                                     temp=data_windowends_states_container{1,i}(n+s*nn-a*nn,j);
                    data_windowends_states_container{1,j}(n+s*nn-a*nn,2:k) = data_backoff_states_container{1,j}(n,2:k)*p_b(no)+...
                        data_windowends_states_container{1,j}(n+s*nn-a*nn,2:k);
%                         zw_disp(data_windowends_states_container{1,j}(n+s*nn-a*nn,2:k));
%                         zw_printf('no: %d, p_a()=%f, p_b()=%f, p_c()=%f, p_d()=%f\n', no, p_a(no),p_b(no),p_c(no),p_d(no));
%                                     zw_printf('12. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+s*nn-a*nn,j, ...
%                                         data_windowends_states_container{1,i}(n+s*nn-a*nn,j),data_backoff_states_container{1,i}(n,j),p_b(no), temp);
                end

                if n+s*nn+1>a*nn
%                         zw_printf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
%                                     temp=data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
                    data_windowends_states_container{1,j}(n+s*nn-a*nn+1,2:k) = data_backoff_states_container{1,j}(n,2:k)*p_c(no)+...
                        data_windowends_states_container{1,j}(n+s*nn-a*nn+1,2:k);
%                                     zw_printf('13. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+s*nn-a*nn+1,j, ...
%                                         data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j),data_backoff_states_container{1,i}(n,j),p_c(no), temp);
                end

%                     else
                if n+s*nn<=(a+1)*nn
                    wait_datawin_ends_temp(n,1) = data_backoff_states_container{1,j}(n,1)*p_a(no)+...
                        wait_datawin_ends_temp(n,1);
                    if n+s*nn>(a)*nn
                        data_windowends_states_container{1,j}(nn-no+1,1) = data_backoff_states_container{1,j}(n,1)*p_d(no)+...
                            data_windowends_states_container{1,j}(nn-no+1,1);
                        continue;
                    else
                        if j~=l+1
                            if j<m+1
                                data_backoff_states_container{1,j+1}(n+s*nn,1:k*2) = data_backoff_states_container{1,j}(n,1)*p_d(no)/k/2+...
                                    data_backoff_states_container{1,j+1}(n+s*nn,1:k*2);
                            else
                                data_backoff_states_container{1,j+1}(n+s*nn,1:k) = data_backoff_states_container{1,j}(n,1)*p_d(no)/k+...
                                    data_backoff_states_container{1,j+1}(n+s*nn,1:k);
                            end
                        end
                    end
                else
                    if n+1*nn>a*nn
%                         zw_printf('case 1, n: %d, j: %d, 1*nn: %d, a*nn: %d\n', n, j, 1*nn, a*nn);
%                             temp=data_windowends_states_container{1,i}(n+1*nn-a*nn,j);
                        data_windowends_states_container{1,j}(n+1*nn-a*nn,1) = data_backoff_states_container{1,j}(n,1)*p_a(no)+...
                            data_windowends_states_container{1,j}(n+1*nn-a*nn,1);
%                             zw_printf('21. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+1*nn-a*nn,j, ...
%                                 data_windowends_states_container{1,i}(n+1*nn-a*nn,j),data_backoff_states_container{1,i}(n,j),p_a(no), temp);
                    else
                        data_windowends_states_container{1,j}(nn-no+1,1) = data_backoff_states_container{1,j}(n,1)*1+...
                            data_windowends_states_container{1,j}(nn-no+1,1);
                        continue;
                    end

                    if n+s*nn>a*nn
%                         zw_printf('case 2, n: %d, j: %d, s*nn: %d, a*nn: %d\n', n, j, s*nn, a*nn);
%                             temp=data_windowends_states_container{1,i}(n+s*nn-a*nn,j);
                        data_windowends_states_container{1,j}(n+s*nn-a*nn,1) = data_backoff_states_container{1,j}(n,1)*p_b(no)+...
                            data_windowends_states_container{1,j}(n+s*nn-a*nn,1);
%                             zw_printf('22. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+s*nn-a*nn,j, ...
%                                 data_windowends_states_container{1,i}(n+s*nn-a*nn,j),data_backoff_states_container{1,i}(n,j),p_b(no), temp);
                    end

                    if n+s*nn+1>a*nn
%                         zw_printf('case 2, n: %d, j: %d, s*nn+1: %d, a*nn: %d\n', n, j, s*nn+1, a*nn);
%                             temp=data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j);
                        data_windowends_states_container{1,j}(n+s*nn-a*nn+1,1) = data_backoff_states_container{1,j}(n,1)*p_c(no)+...
                            data_windowends_states_container{1,j}(n+s*nn-a*nn+1,1);
%                             zw_printf('23. n: %d, winends_cont{1,%d}(%d,%d) %g=%g*%g+%g\n', n, i, n+s*nn-a*nn+1,j, ...
%                                 data_windowends_states_container{1,i}(n+s*nn-a*nn+1,j),data_backoff_states_container{1,i}(n,j),p_c(no), temp);
                    end
                end

            end
%                     end

%              toc(tt);

%                 end
%                 toc(t3);
%             end
%                 zw_printf('iter: %d count: %d\n',j, size(row,1));

       end   

       data_p_save_container{zw} = p_save;

       toc(t2);

%            return;


%====================================================================================


%                toc;
%             zw_printf('data_backoff_states_container{1,%d} after...\n', i);
%             data_backoff_states_container{1,i}
%             for tt=k:-1:k-2
% tic;
% 
%                 size(find(data_backoff_states_container{1,i}(:,tt)))
%                 [row_temp]=find(data_backoff_states_container{1,i}(:,tt));
%                 
%                 zw_printf('%d\n',row_temp);
%                 for zw=row_temp
%                     zw_printf('%g\n', data_backoff_states_container{1,i}(zw,tt));
%                 end
%                 
%                 
%     toc;
%             end
%         continue;

        if SYMBOL_VERIFICATION
            zw_printf('data_backoff_states_container{1,%d} after...\n', i);
            zw_disp(simplify(data_backoff_states_container{1,i}));

            zw_printf('data_windowends_states_container{1,%d} after...\n', i);
            zw_disp(simplify(data_windowends_states_container{1,i}));

            zw_disp('sum_data_backoff_states:');
            zw_disp(simplify(sum_data_backoff_states));

            wait_datawin_ends_temp = sym(zeros((a-s+1)*nn,1)); %reset wait_datawin_ends_temp
            wait_datawin_ends_temp2 = sym(zeros((a-s+1)*nn,1)); %reset wait_datawin_ends_temp
        else
%                 zw_printf('data_backoff_states_container{1,%d} after...\n', i);
%                 zw_disp((data_backoff_states_container{1,i}));
% 
%                 zw_printf('data_windowends_states_container{1,%d} after...\n', i);
%                 zw_disp((data_windowends_states_container{1,i}));
% 
            zw_disp('sum_data_backoff_states:');
            zw_disp((sum_data_backoff_states));
% 
            wait_datawin_ends_temp2 = zeros((a-s+1)*nn,1); %reset wait_datawin_ends_temp
        end

%             return;
%             tttt=tic;

        throughput_container_temp{zw}=sum(wait_datawin_ends_temp)*s;
        
        data_energy_container_temp{zw}=0;
        
        for n=1:nn:(a-s+1)*nn
            data_energy_container_temp{zw} = sum(wait_datawin_ends_temp(n:n+nn-1)*(ceil(n/nn)+s+w_a)) + data_energy_container_temp{zw};
        end
        
        fprintf('%g\n', data_energy_container_temp{zw});
        

        offset = s*nn;
        for n=1:(a-s)*nn
        %             zw_disp(m);
        %             zw_disp(n);
%                 if 1%data_backoff_states_container{1,i}(n,1)~=0

            if mod(n,nn)==0
                no = 1;
            else
                no = nn - mod(n,nn) + 1;
            end

%                 zw_printf('n: %d\n', n);

            if wait_datawin_ends_temp(n,1)~=0

                if NEWMC == 3
                    if no~=1
                        p_a(no)=p_save(n+offset+1,1);
                        p_b(no)=p_save(n+offset+1,2);
                        p_c(no)=p_save(n+offset+1,3);

                        if p_a(no)+p_b(no)+p_c(no)~=1
                            p_a(no)=1;
                            p_b(no)=0;
                            p_c(no)=0;                                
                        end

                    else
                        p_a(no)=1;
                        p_b(no)=0;
                        p_c(no)=0;
                    end
                end


                if n+1*nn<=(a-s+1)*nn
%                         zw_printf('case 1: %d to %d with %g\n', n+offset, n+offset+1*nn, p_a(no));
                    wait_datawin_ends_temp(n+1*nn,1) = wait_datawin_ends_temp(n,1)*p_a(no)+...
                        wait_datawin_ends_temp(n+1*nn,1);
                end 

                if n+s*nn<=(a-s+1)*nn
%                         zw_printf('case 2: %d to %d with %g\n', n+offset, n+offset+s*nn, p_b(no));
                    wait_datawin_ends_temp(n+s*nn,1) = wait_datawin_ends_temp(n,1)*p_b(no)+...
                        wait_datawin_ends_temp(n+s*nn,1);
                    wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp2(n,1) + ...
                        wait_datawin_ends_temp(n,1)*p_b(no)*(s-1);
%                         zw_printf('case 2: %g\n', wait_datawin_ends_temp2(n,1));
                else
%                         zw_printf('case 4: %d to %d with %g\n', n+offset, (a-s+1)*nn-no+1+offset, p_b(no));
                    wait_datawin_ends_temp((a-s+1)*nn-no+1,1) = wait_datawin_ends_temp(n,1)*p_b(no)+...
                        wait_datawin_ends_temp((a-s+1)*nn-no+1,1);

                    wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp2(n,1) + ...
                        wait_datawin_ends_temp(n,1)*p_b(no)*(((a-s+1)*nn-no+1-n)/nn-1);
%                         zw_printf('case 4: n: %d, no: %d, %d, %g\n', n, no, ((a-s+1)*nn-no+1-n)/nn-1, wait_datawin_ends_temp2(n,1));
                end

                if n+s*nn+1<=(a-s+1)*nn
%                         zw_printf('case 3: %d to %d with %g\n', n+offset, n+offset+s*nn+1, p_c(no));
                    wait_datawin_ends_temp(n+s*nn+1,1) = wait_datawin_ends_temp(n,1)*p_c(no)+...
                        wait_datawin_ends_temp(n+s*nn+1,1);
                    wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp2(n,1) + ...
                        wait_datawin_ends_temp(n,1)*p_c(no)*(s-1);
%                         zw_printf('case 3: n: %d %g\n', n, wait_datawin_ends_temp2(n,1));
                else
%                         zw_printf('case 5: %d to %d with %g\n', n+offset, (a-s+1)*nn-no+1+offset, p_c(no));
                    wait_datawin_ends_temp((a-s+1)*nn-no+1,1) = wait_datawin_ends_temp(n,1)*p_c(no)+...
                        wait_datawin_ends_temp((a-s+1)*nn-no+1,1);
                    wait_datawin_ends_temp2(n,1) = wait_datawin_ends_temp2(n,1) + ...
                        wait_datawin_ends_temp(n,1)*p_c(no)*(((a-s+1)*nn-no+1-n)/nn-1);
%                         zw_printf('case 5: n: %d, no: %d, %d, %g\n', n, no, ((a-s+1)*nn-no+1-n)/nn-1, wait_datawin_ends_temp2(n,1));

                end

            end

        end


%             return;

%             toc(tttt);


         if SYMBOL_VERIFICATION

%             zw_disp('wait_datawin_ends_temp:');
%             zw_disp(simplify(wait_datawin_ends_temp));
         else
%                 zw_disp('wait_datawin_ends_temp:');
%                 zw_disp((wait_datawin_ends_temp));
         end


        wait_datawin_ends = wait_datawin_ends + wait_datawin_ends_temp;
        wait_datawin_ends2 = wait_datawin_ends2 + wait_datawin_ends_temp2;

%             zw_disp('data_backoff_states_container:');
%             zw_disp((data_backoff_states_container{1,i}));
%             zw_disp('data_windowends_states_container:');
%             zw_disp((data_windowends_states_container{1,i}));
%         
%         
        if NEWMC==3
            for j=1:l+1
                sum_data_windowends_states = sum_data_windowends_states + sum(sum(data_windowends_states_container{1,j}));
            end
            zw_disp('sum_data_windowends_states:');
            zw_disp(sum_data_windowends_states);

            break;
        end

        sum_data_windowends_states = sum_data_windowends_states + sum(sum(data_windowends_states_container{1,i}));
        zw_disp('sum_data_windowends_states:');
        zw_disp(sum_data_windowends_states);
%         
%         zw_disp('wait_datawin_ends:');
%         zw_disp(simplify(wait_datawin_ends));

    end


    toc(tt);


    zw_printf('=======DATA summary for active nodes: %d=======\n', nn-zw+1);

    %compute the p2f here.

    Ct=data_backoff_states_container{1,l+1}(:,1);
%                 zw_disp((data_backoff_states_container{1,i}));
    [row]=find(Ct>PRECISION);


    for n=row.'
        if mod(n,nn)==0
            no = 1;
        else
            no = nn - mod(n,nn) + 1;
        end

        if NEWMC == 3
            p_d(no)=p_save(n,4);
        end

        if n+s*nn<=a*nn
            p2f(nn+1-zw,nn-no+1-zw+1) = p2f(nn+1-zw,nn-no+1-zw+1) + Ct(n,1)*p_d(no);
%                 zw_printf('%d with %g\n', n, p_d(no));
%                 zw_printf('nn-no+1-zw+1: %d, p2f(nn+1-zw,nn-no+1-zw+1): %g\n', nn-no+1-zw+1, p2f(nn+1-zw,nn-no+1-zw+1));
        end
    end

%     zw_printf('p2f:\n');
%     zw_disp(p2f);


%     cc=0;
%     for n=nn+1-zw:-1:1
%         ccc=0;
%         for k=1:nn-cc
% %                     zw_printf('%d, %d\n', n-1,cc+k-1);
% %                     zw_printf('%d, %d\n', nn-cc-1,ccc);
%             pex(n,cc+k)=nchoosek(nn-cc-1,ccc)*P^(nn-cc-1-ccc)*(1-P)^(ccc);
% %                 zw_printf('%d to %d, %d*P^%d*(1-P)^%d\n', n-1, cc+k, nchoosek(nn-cc-1,ccc), nn-cc-1-ccc, ccc);
%             ccc=ccc+1;
%         end
%         cc=cc+1;
%     end
%     zw_printf('pex:\n');
%     zw_disp(pex);



%     return;


    if SYMBOL_VERIFICATION
        zw_disp('wait_datawin_ends:');
        zw_disp(simplify(wait_datawin_ends));

        zw_disp('wait_datawin_ends2:');
        zw_disp(simplify(wait_datawin_ends2));

        zw_disp('sum_data_backoff_states:');
        sum_data_backoff_states = simplify(sum_data_backoff_states);
        zw_disp(sum_data_backoff_states);   

        zw_disp('sum_data_windowends_states:');
        sum_data_windowends_states = simplify(sum_data_windowends_states);
        zw_disp(sum_data_windowends_states);

        sum_wait_datawin_ends = simplify(sum(wait_datawin_ends)+ sum(wait_datawin_ends2));
        zw_disp('sum_wait_datawin_ends:');
        zw_disp(sum_wait_datawin_ends);
    else

%             zw_disp('wait_datawin_ends:');
%             zw_disp(wait_datawin_ends);
%             zw_disp((wait_datawin_ends((a-s)*nn+1:(a-s+1)*nn,1)));
        p2s(nn-zw+1,:) = circshift(wait_datawin_ends((a-s)*nn+1:(a-s+1)*nn,1),nn+1-zw)';

%             cc=0;
%             for n=nn+1-zw:-1:1
%                 ccc=0;
%                 for k=1:nn-cc
%     %                     zw_printf('%d, %d\n', n-1,cc+k-1);
%     %                     zw_printf('%d, %d\n', nn-cc-1,ccc);
%                     pnx(n,cc+k)=(1-P)*nchoosek(nn-cc-1,ccc)*P^(nn-cc-1-ccc)*(1-P)^(ccc);
%                         zw_printf('%d to %d, (1-P)*%d*P^%d*(1-P)^%d\n', n, cc+k, nchoosek(nn-cc-1,ccc), nn-cc-1-ccc, ccc);
%                     ccc=ccc+1;
%                 end
%                 cc=cc+1;
%             end

% %         pnx = pex*(1-P);
        
%             zw_printf('pnx:\n');
%             zw_disp(pnx);            

%             data_success_probs = data_success_probs + circshift(wait_datawin_ends((a-s)*nn+1:(a-s+1)*nn,1),nn+1-zw);
%             zw_disp(circshift(wait_datawin_ends((a-s)*nn+1:(a-s+1)*nn,1),nn+1-zw));

%             zw_disp('wait_datawin_ends2:');%(a-s+1)*nn
%             zw_disp((wait_datawin_ends2()));
% 
%             zw_disp('sum_data_backoff_states:');
%             zw_disp(sum_data_backoff_states);   
% 
%             zw_disp('sum_data_windowends_states:');
%             zw_disp(sum_data_windowends_states);

        sum_wait_datawin_ends = sum_wait_datawin_ends + sum(wait_datawin_ends) + sum(wait_datawin_ends2);
        zw_disp('sum_wait_datawin_ends:');
        zw_disp(sum_wait_datawin_ends);

        data_windowends_probs_temp = zeros(nn*s+1,1);

        for i=1:l+1

            data_windowends_probs_temp = data_windowends_probs_temp + sum(data_windowends_states_container{1,i}')';

        end

        j = size(data_windowends_probs_temp,1);
        data_windowends_probs_temp(j,:)=[];
%             zw_disp('data_windowends_probs_temp');
%             zw_disp(data_windowends_probs_temp);
        data_windowends_probs = zeros(nn,1);
        for i=1:nn:j-1
             data_windowends_probs = data_windowends_probs + data_windowends_probs_temp(i:i+nn-1,1);
        end

        p2e(nn-zw+1,:) = circshift(data_windowends_probs,nn+1-zw)';

%         zw0 = 0;
%         zw1 = 0;
%         for i=1:l+1
%             
%             zw0 = zw0 + (sum(data_backoff_states_container{1,i}(:,1)));
%             zw1 = zw1 + (sum(sum(data_backoff_states_container{1,i})));
%             
%             
%         end
%             tau = zw0/zw1;
%             zw_printf('tau: %g\n', tau);

    end

%     toc(t);
    sum_states = sum_data_backoff_states + sum_data_windowends_states + sum_wait_datawin_ends;
%         zw_printf('all: %g\n', all);

%     zw_disp('data_success_probs:');
%     zw_disp(data_success_probs);
%     
%     
%     zw_disp('data_windowends_probs:');
%     zw_disp(data_windowends_probs);




%         zw_printf('sum p2e:\n');
%         zw_disp(sum(p2e(4,:)));
%         zw_printf('sum p2f:\n');
%         zw_disp(sum(p2f(4,:)));
%         zw_printf('sum p2s:\n');
%         zw_disp(sum(p2s(4,:)));


    for i=1:nn
        haha=sum(p2s(i,:))+sum(p2e(i,:))+sum(p2f(i,:));
        zw_printf('%g\n', haha);
    end

%     pwx = pex*P;
% %         zw_printf('pwx:\n');
% %         zw_disp(pwx);
% 
%     p1 = p1 + p2e*pex;
%     p2 = p2 + p2f*pex;
%     p3 = p3 + p2s*pnx;
%     p4 = p4 + p2s*pwx;
    data_energy_container_temp{zw} = sum(sum(p2e))*(w_d+w_a) + data_energy_container_temp{zw};
    data_energy_container_temp{zw} = sum(sum(p2f))*(w_d+w_a) + data_energy_container_temp{zw};
    data_pe_container{zw} = p2e;
    data_pf_container{zw} = p2f;
    data_ps_container{zw} = p2s;
    data_sumstate_container{zw} = sum_states;
    
end

