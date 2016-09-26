function [  ] = build_waiting_states( window_type, zw )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    global nn;
    global w_d;    
    global w_a;    
    global s_d;    
    global s_a;    
    global data_p_save_container;
    global atim_p_save_container;
    global WINDOW_ATIM;
    global WINDOW_DATA;
    global atim_wait_container;
    global data_wait_container;
    
   
    if window_type == WINDOW_DATA
        a = w_d;
        s = s_d;
        p_save_container = data_p_save_container;
    else
        a = w_a;
        s = s_a;
        p_save_container = atim_p_save_container;
    end


    wait_temp = zeros((a+1)*nn,1); 
    wait_temp2 = zeros((a+1)*nn,1);

    wait_temp(zw) = 1;

    for n=1:(a)*nn
    %             zw_disp(m);
    %             zw_disp(n);
%                 if 1%atim_backoff_states_container{1,i}(n,1)~=0

        if mod(n,nn)==0
            no = 1;
        else
            no = nn - mod(n,nn) + 1;
        end

%                 zw_printf('no: %d\n', no);

        if wait_temp(n,1)~=0

            if no~=1
                p_a(no)=p_save_container{1,zw+1}(n+1,1);
                p_b(no)=p_save_container{1,zw+1}(n+1,2);
                p_c(no)=p_save_container{1,zw+1}(n+1,3);

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

            if n+1*nn<=(a+1)*nn
                wait_temp(n+1*nn,1) = wait_temp(n,1)*p_a(no)+...
                    wait_temp(n+1*nn,1);
%                         zw_printf('case 1: %d to %d with %g\n', n, n+1*nn, p_a(no));
            end 

            if n+s*nn<=(a+1)*nn
                wait_temp(n+s*nn,1) = wait_temp(n,1)*p_b(no)+...
                    wait_temp(n+s*nn,1);
%                         zw_printf('case 2: %d to %d with %g\n', n, n+s*nn, p_b(no));
                wait_temp2(n,1) = wait_temp2(n,1) + ...
                    wait_temp(n,1)*p_b(no)*(s-1);
            else
                wait_temp((a+1)*nn-no+1,1) = wait_temp(n,1)*p_b(no)+...
                    wait_temp((a+1)*nn-no+1,1);
%                         zw_printf('case 4: %d to %d with %g\n', n, (a+1)*nn-no+1, p_b(no));
                wait_temp2(n,1) = wait_temp2(n,1) + ...
                    wait_temp(n,1)*p_b(no)*(((a+1)*nn-no+1-n)/nn-1);
            end

            if n+s*nn+1<=(a+1)*nn
                wait_temp(n+s*nn+1,1) = wait_temp(n,1)*p_c(no)+...
                    wait_temp(n+s*nn+1,1);
%                         zw_printf('case 3: %d to %d with %g\n', n, n+s*nn+1, p_c(no));
                wait_temp2(n,1) = wait_temp2(n,1) + ...
                    wait_temp(n,1)*p_c(no)*(s-1);
            else
                wait_temp((a+1)*nn-no+1,1) = wait_temp(n,1)*p_c(no)+...
                    wait_temp((a+1)*nn-no+1,1);
%                         zw_printf('case 5: %d to %d with %g\n', n, (a+1)*nn-no+1, p_c(no));
                wait_temp2(n,1) = wait_temp2(n,1) + ...
                    wait_temp(n,1)*p_c(no)*(((a+1)*nn-no+1-n)/nn-1);

            end

        end

    end

    
%     zw_printf('wait_temp:\n');
%     zw_disp(wait_temp);
%     zw_printf('wait_temp2:\n');
%     zw_disp(wait_temp2);


    pw = zeros(nn,nn);

    pw(nn+1-zw,1:nn+1-zw)=wait_temp((a)*nn+zw:(a+1)*nn,1)';
    
    
%             zw_printf('pw:\n');
%             zw_disp(pw);

%     pws = zeros(nn,nn);
% 
%     cc=0;
%     for n=nn+1-zw:-1:1
%         ccc=0;
%         for k=1:nn-cc
% %                     zw_printf('%d, %d\n', n-1,cc+k-1);
% %                     zw_printf('%d, %d\n', nn-cc-1,ccc);
%             pws(n,cc+k)=P*nchoosek(nn-cc-1,ccc)*P^(nn-cc-1-ccc)*(1-P)^(ccc);
% %                     zw_printf('%d to %d, P*%d*P^%d*(1-P)^%d\n', n-1, cc+k-1, nchoosek(nn-cc-1,ccc), nn-cc-1-ccc, ccc);
%             ccc=ccc+1;
%         end
%         cc=cc+1;
%     end
% %             zw_printf('pws:\n');
% %             zw_disp(pws);
% 
%     p5 = p5 + pw*pws;
% %             zw_printf('t_matrix:\n');
% %             zw_disp(t_matrix);

    sum_states = sum(wait_temp) + sum(wait_temp2);
    
    if window_type == WINDOW_DATA
        data_wait_container{1,zw} = pw;
        data_wait_container{2,zw} = sum_states;
    else
        atim_wait_container{1,zw} = pw;
        atim_wait_container{2,zw} = sum_states;
    end

    
end

