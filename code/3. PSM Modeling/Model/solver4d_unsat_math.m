
function f=solver4d_unsat_math(xvect)
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


    global PRECISION_INIT;
    global SYMBOL_VERIFICATION;
    global DEBUG;
    global TEST;
    global WINDOW_ATIM;
    global WINDOW_DATA;
    global lambda
    
    global p125;
    global Ls;
    global mu;
    
    if TEST==2
        global P;
    end

    ENABLE_TICTOC_DETAILED = 0;
    ENABLE_TICTOC = 1;



%     for i=1:l_a+1 %generate all A's (i.e. A1 A2 ...) and B's
%         %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
%         %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
%         eval(['global A',num2str(i),';']);
%         eval(['global B',num2str(i),';']);
%     %     eval(['A',num2str(i),'=zeros(1,W*2^(i-1))']);
%     %     B{i}=eval(['A',num2str(i)]);
%     end
%     for i=1:l_d+1 %generate all A's (i.e. A1 A2 ...) and B's
%         %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
%         %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
%         eval(['global C',num2str(i),';']);
%         eval(['global D',num2str(i),';']);
%     end

    for i=1:nn
        eval(['init', num2str(i), '=xvect(', num2str(i), ');']);
    end
    if TEST==0
        eval(['P=xvect(', num2str(i+1), ');']);
    end

    
    fprintf('New iteration...\n');
%      t=tic;
    p1 = zeros(nn,nn); 
    p2 = zeros(nn,nn); 
    p3 = zeros(nn,nn);
    p4 = zeros(nn,nn);
    p5 = zeros(nn,nn);
    p6 = zeros(nn,nn);
        
    all = 0;
      
    for zw=1:nn%nn%
        zw_printf('=======Begin with ATIM for active nodes: %d=======\n', nn-zw+1);
        
%         atim_windowends_probs = zeros(nn,1); 
        %store the window ends probability for different number of nodes, from the first to the last line, are the prob. of one to nn nodes
        %not transmit successfully, i.e. the prob. of one to nn nodes left due to window ends.
        init = eval(['init',num2str(zw),';']);
%         disp(init);
        [pe, pf, ps, sum_temp, ~] = handle_atimwin( zw, init );

        all = all + sum_temp;
        
        zw_printf('pe:\n');
        zw_disp(pe);
        zw_printf('pf:\n');
        zw_disp(pf);
        zw_printf('ps:\n');
        zw_disp(ps);
        zw_printf('sum_temp:\n');
        zw_disp(sum_temp);
        
        pef = pe + pf;
        
%         return;
        
        %get p2s, p2e and p2f
        p2e = zeros(nn,nn);
        p2f = zeros(nn,nn);
        p2s = zeros(nn,nn);
        
        zw_printf('====Now begin with DATA window====\n');
        
        for zw2=1:nn
%             [p2e_temp, p2f_temp, p2s_temp, sum_temp] = handle_datawin( zw2, 0.21 );
            [p2e_temp, p2f_temp, p2s_temp, sum_temp, ~, ~] = handle_datawin( zw2, ps(nn+1-zw, nn+1-zw2) );
            p2e = p2e + p2e_temp;
            p2f = p2f + p2f_temp;
            p2s = p2s + p2s_temp;
            all = all + sum_temp;
    
%         zw_printf('p2e_temp:\n');
%         zw_disp(p2e_temp);
%         zw_printf('p2f_temp:\n');
%         zw_disp(p2f_temp);
%         zw_printf('p2s_temp:\n');
%         zw_disp(p2s_temp);

% return;
        end
        
%         continue;
        
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(p2e);
        p2e = temp_matrix;
        
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(p2f);
        p2f = temp_matrix;
        
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(p2s);
        p2s = temp_matrix;
        
        pxn = compute_pxn(zw, P);
        
        p1 = p1 + (p2e)*pxn;
        p2 = p2 + (p2f)*pxn;
        p3 = p3 + (p2s)*(1-P)*pxn;
        p4 = p4 + (p2s)*(P)*pxn;
        
%         zw_printf('p1:\n');
%         zw_disp(p1);
%         zw_printf('p2:\n');
%         zw_disp(p2);
%         zw_printf('p3:\n');
%         zw_disp(p3);
%         zw_printf('p4:\n');
%         zw_disp(p4);
%         zw_printf('all: %g\n', all);
        
        
        
        
        %get pew
        pew = zeros(nn,nn);
%         DEBUG = 1;
        for zw2=1:nn
            [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pef(nn+1-zw, nn+1-zw2) );
%             [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, 1 );
            pew = pew + pew_temp;
            all = all + sum_temp;
    
%         zw_printf('pew_temp:\n');
%         zw_disp(pew_temp);
%         zw_printf('sum_temp:\n');
%         zw_disp(sum_temp);

% return;

        end
        
        
        
% return;
        
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(pew);
        pew = temp_matrix;
        
%         zw_printf('pew:\n');
%         zw_disp(pew);

        
        p5 = p5 + (pew)*pxn;
        
%         return;
        
%         zw_printf('p2e:\n');
%         zw_disp(sum(p2e));
%         zw_printf('p2f:\n');
%         zw_disp(sum(p2f));
%         zw_printf('p2s:\n');
%         zw_disp(sum(p2s));
%         
%         return;
        
%         zw_printf('p2e:\n');
%         zw_disp(ps*p2e);
%         zw_printf('p2f:\n');
%         zw_disp(ps*p2f);
%         zw_printf('p2s:\n');
%         zw_disp(ps*p2s);        
        
        

    
    end
%     return;
%     toc(t);
    
    DEBUG=0;
    zw_printf('p1:\n');
    zw_disp(p1);
    zw_printf('p2:\n');
    zw_disp(p2);
    zw_printf('p3:\n');
    zw_disp(p3);
    zw_printf('p4:\n');
    zw_disp(p4);
    zw_printf('p5:\n');
    zw_disp(p5);
    DEBUG=0;
    
   
%     ones(1,nn)*p1
%     ones(1,nn)*p2
%     ones(1,nn)*p3

    b=sum(p4);
    zw_printf('Input for waiting queue...\n');
    zw_disp(b');
    
%     return;
    
%     continue;
    
    zw_printf('Handling waiting queue...\n');


    pw = zeros(nn,nn);
    for zw=1:nn
%             [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pef(nn+1-zw, nn+1-zw2) );
        [pw_temp, sum_temp] = waiting_states( WINDOW_ATIM, zw, 1 );
        pw = pw_temp;
%             all = all + sum_temp;

%         zw_printf('pw_temp:\n');
%         zw_disp(pw_temp);
%         zw_printf('sum_temp:\n');
%         zw_disp(sum_temp);

        p2w = zeros(nn,nn);
        for zw2=1:nn
            [p2w_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pw(nn+1-zw, nn+1-zw2) );
%             [p2w_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, 1 );
            p2w = p2w + p2w_temp;
%                 all = all + sum_temp;

        zw_printf('p2w_temp:\n');
        zw_disp(p2w_temp);
%         zw_printf('sum_temp:\n');
%         zw_disp(sum_temp);


        end
        
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(p2w);
        p2w = temp_matrix;

        pxn = compute_pxn(zw, P);
        p6 = p6 + p2w*P*pxn;

    end 
        
%         temp_matrix = zeros(nn,nn);
%         temp_matrix(nn+1-zw,:) = sum(pw);
%         pew = temp_matrix;
        
       
        
    A=p6'-eye(nn);
%         b=[-1;-1;-1;-1];%initial input prob. (from successful tx of ATIM/DATA part) for different # active nodes, from 0 to nn-1.
    b=-1*b';
    selfloop_final_inputprobs = A\b; %final input prob. (from successful tx of ATIM/DATA part) for different # active nodes, from 0 to nn-1
    DEBUG=0;
    zw_printf('selfloop_final_inputprobs:\n');
    zw_disp(selfloop_final_inputprobs);
    DEBUG=0;

    
    p6 = zeros(nn,nn);
    pw = zeros(nn,nn);
    for zw=1:nn
%             [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pef(nn+1-zw, nn+1-zw2) );
        [pw_temp, sum_temp] = waiting_states( WINDOW_ATIM, zw, selfloop_final_inputprobs(nn+1-zw) );
        pw = pw_temp;
        all = all + sum_temp;

%         zw_printf('pew_temp:\n');
%         zw_disp(pew_temp);
%         zw_printf('sum_temp:\n');
%         zw_disp(sum_temp);

        p2w = zeros(nn,nn);
        for zw2=1:nn
            [p2w_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pw(nn+1-zw, nn+1-zw2) );
%             [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, 1 );
            p2w = p2w + p2w_temp;
            all = all + sum_temp;

%         zw_printf('p2w_temp:\n');
%         zw_disp(p2w_temp);
%         zw_printf('sum_temp:\n');
%         zw_disp(sum_temp);

        end
        temp_matrix = zeros(nn,nn);
        temp_matrix(nn+1-zw,:) = sum(p2w);
        p2w = temp_matrix;

        pxn = compute_pxn(zw, P);
        p6 = p6 + p2w*(1-P)*pxn;

    end 

    zw_printf('p6:\n');
    zw_disp(p6);
    DEBUG=1;
    zw_printf('all: %g, P: %g, ', all, P);
    DEBUG=0;
%         return;
%     end
   

    init_new = sum(p1+p2+p3+p5+p6);
%     fprintf('init_new:\n');
%     disp(init_new');
    if TEST==0

        p125=compute_p125(P);
        DEBUG=0;
        zw_printf('p125:\n');
        zw_disp(p125);
        DEBUG=0;
        L=(w_a+w_d)/100;
        Ls=(p125-eye(nn))\(-L*ones(nn,1));

%         mu=sum(p3+p6)/sum(sum(p3+p6))*p125*Ls+L;
        mu=init_new/sum(init_new)*p125*Ls+L;
%         mu=sum(p1+p2+p5)/sum(sum(p1+p2+p5))*p125*Ls+L;
        DEBUG=1;
        zw_printf('mu: %f, ', 1000/mu);
        new_P = 1-lambda/(1000/mu);%???? What if using M/G/1/K queue?
        zw_printf('new_P: %f\n', new_P);
        DEBUG=0;
    else
        DEBUG=1;
        zw_printf('\n');
        DEBUG=0;
        
    end


    f(1) = 1 - all;
    for i=1:nn
        tmp = sprintf('f(%d+1) = init_new(%d)-init%d;', i, i, nn+1-i);
        eval(tmp);
    end    
    if TEST==0
        f(nn+2) = new_P - P;
    end
%     f(2) = computed_tauas(1) - taua1;
%     f(3) = computed_tauas(2) - taua2;
%     f(4) = computed_tauas(3) - taua3;
%     f(5) = computed_tauas(4) - taua4;
%     f(6) = init_new(1) - init4;
%     f(7) = init_new(2) - init3;
%     f(8) = init_new(3) - init2;
%     f(9) = init_new(4) - init1;
    
    
    












