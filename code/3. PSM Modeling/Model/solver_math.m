
function f=solver_math(xvect)
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

global atim_backoff_states_container;% atim_backoff_states_container is a cell for containing all A's (i.e. A1 A2 ...)
global atim_windowends_states_container;% atim_windowends_states_container is a cell for containing all B's (i.e. B1 B2 ...)
global sum_atim_backoff_states;% sum_atim_backoff_states is for storing the sum of atim_backoff_states_container
global sum_atim_windowends_states;% sum_atim_backoff_states is for storing the sum of atim_windowends_states_container
global wait_atimwin_ends; %waiting for the end of the ATIM window
global wait_atimwin_ends_temp;
global sum_wait_atimwin_ends; %for storing the sum of wait_atimwin_ends


for i=1:l+1 %generate all A's (i.e. A1 A2 ...) and B's
    %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    eval(['global A',num2str(i),';']);
    eval(['global B',num2str(i),';']);
%     eval(['A',num2str(i),'=zeros(1,W*2^(i-1))']);
%     B{i}=eval(['A',num2str(i)]);
end
% syms p;
% syms bi0;
% 
% 
% b00 = symsum(bi0,0,m);
% W0 = W;

    p = xvect(1);
    init = xvect(2);
%     tau = xvect(3);

%the followings are for the Bianchi's model, using method in constructive
%manner.


    for i=1:l+1
        atim_backoff_states_container{i}=eval(['A',num2str(i)]);
        atim_windowends_states_container{i}=eval(['B',num2str(i)]);
    end
    
%     p=0.2;
    wait_atimwin_ends = 0;
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
                            if n+s*nn<=a*nn
                                atim_backoff_states_container{1,i}(n+s*nn,j) = Ct(n,1)*p/k+atim_backoff_states_container{1,i}(n+s*nn,j);
                            end
                        end
                    end
                end            
            end
            
%             
%             disp('atim_backoff_states_container{1,i}');
%             disp(atim_backoff_states_container{1,i});
            
           for j=k:-1:1
                
                for n=1:a*nn
                %             disp(m);
                %             disp(n);
                    if atim_backoff_states_container{1,i}(n,j)~=0
                        
                        if j~=1
                            if n<=a*nn-s+1
                                sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1)*p);
                                %(s-1)*p represents the freezing states for s-1 time slots due to the transmission of some other nodes
                            else
                                sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a*nn-n)*p);
                                %here since if ATIM window ends, the freezing time may not be s-1 long.
                            end
                        else %j==1
                            if n<=a*nn-s+1
                                sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(s-1));
                            else
                                sum_atim_backoff_states = sum_atim_backoff_states + atim_backoff_states_container{1,i}(n,j)*(1+(a*nn-n)*1);
                            end
                        end
%                         disp('sum_atim_backoff_states:');
%                         disp(simplify(sum_atim_backoff_states));
                        
                        if j>1 %except the first column
                %                     disp('m>1');
                            if n+1*nn<=a*nn
                                atim_backoff_states_container{1,i}(n+1*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*(1-p)+atim_backoff_states_container{1,i}(n+1*nn,j-1);
                            end          
                            if n+s*nn<=a*nn
                                atim_backoff_states_container{1,i}(n+s*nn,j-1) = atim_backoff_states_container{1,i}(n,j)*p+atim_backoff_states_container{1,i}(n+s*nn,j-1);
                            end
                        end

                        if n+1*nn>a*nn
                            atim_windowends_states_container{1,i}(n+1*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*(1-p)+atim_windowends_states_container{1,i}(n+1*nn-a*nn,j);
            %                         disp(Bt(m-1,n+1-k));
                        end

                        if n+s*nn>a*nn
                            atim_windowends_states_container{1,i}(n+s*nn-a*nn,j) = atim_backoff_states_container{1,i}(n,j)*p+atim_windowends_states_container{1,i}(n+s*nn-a*nn,j);
            %                         disp(Bt(m-1,n+s-k));
                        end
                    end
                end
           end
           

        wait_atimwin_ends_temp = zeros(a-s+1,1); %reset wait_atimwin_ends_temp
        
         for n=1:a*nn-s+1
        %             disp(m);
        %             disp(n);
            if 1%atim_backoff_states_container{1,i}(n,1)~=0
                if n==1
                    wait_atimwin_ends_temp(n,1) = atim_backoff_states_container{1,i}(n,1)*(1-p);
                else
                    wait_atimwin_ends_temp(n,1) = wait_atimwin_ends_temp(n-1,1)+ atim_backoff_states_container{1,i}(n,1)*(1-p);
                end
            end
         end
         
%         disp('wait_atimwin_ends_temp:');
%         disp(simplify(wait_atimwin_ends_temp));

        
        wait_atimwin_ends = wait_atimwin_ends + wait_atimwin_ends_temp;
        
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

%         disp('sum_atim_backoff_states:');
% %         sum_atim_backoff_states = simplify(sum_atim_backoff_states);
%         disp(sum_atim_backoff_states);   
%         
%         disp('sum_atim_windowends_states:');
% %         sum_atim_windowends_states = simplify(sum_atim_windowends_states);
%         disp(sum_atim_windowends_states);
%         
% %         sum_wait_atimwin_ends = simplify(sum(wait_atimwin_ends));
        sum_wait_atimwin_ends = sum(wait_atimwin_ends);
%         disp('sum_wait_atimwin_ends:');
%         disp(sum_wait_atimwin_ends);
        
        all = sum_atim_backoff_states + sum_atim_windowends_states + sum_wait_atimwin_ends;
        fprintf('all: %g\n', all);
        
        zw0 = 0;
        zw1 = 0;
        for i=1:l+1
            
            zw0 = zw0 + (sum(atim_backoff_states_container{1,i}(:,1)));
            zw1 = zw1 + (sum(sum(atim_backoff_states_container{1,i})));
            
            
        end
        tau = zw0/zw1;
        fprintf('tau: %g, init: %g, p: %g\n', tau, init, p);
     
    
    f(1) =  1 - (1-tau)^(N-1) - p;

    f(2) = 1/all - init;
    












