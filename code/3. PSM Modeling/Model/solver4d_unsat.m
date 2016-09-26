% m=4;
% n=3;
function [] = solver4d_unsat( n, payload, tr, atim, data, p0 )
% clear all;
% clear all;
% clear all;
% 
% if 1
%     
% % n = 4;
% % tr = 40;
% % payload = 15;
% % atim = 10;
% % data = 100 - atim;
% 
% n = 30;
% tr = 40;
% payload = 1500;
% atim = 100;
% data = 2000 - atim;
% 
% p0 = 0;
% 
% 
    global TEST;
%     TEST = 2;

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
    zw_time(1:6)=0;
    zw_time=clock();
    outname=['run-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
        num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
        '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];
    % diary(outname);
    % diary on;

    zw_disp('**************************Markov Chain Model***************************');
    zw_disp('**************************parameters***************************');
    global PRECISION_INIT;
    PRECISION_INIT=0;%1e-8;%
    % fprintf('PRECISION_INIT: %f\n', PRECISION_INIT);
    global DEBUG;
    DEBUG=0;
    global tau;
    global tauc;

    global lambda
    lambda=tr;

    if TEST==2
        global P;
    end
    P=p0;
    fprintf('P: %f\n', P);

    global p125;
    global Ls;
    global mu;

    %for ATIM window
    %=========================================
    Lo=0;%20;

    PAYLOAD = payload;%bytes
    HEDEAR = 40;%bytes, mac header, tcp/udp header and ip header
    PREAMBLE = 4;%40 us
    LDATA=PREAMBLE+ceil(HEDEAR*8/54/10)+ceil(PAYLOAD*8/54/10);%ceil(846*8/54/10)
    ATIM=30;%bytes
    LATIM=PREAMBLE+ceil(ATIM*8/54/10);

    W_a=16; %contention window size 32
    m_a=5; %the CW max. 5
    l_a=5; %the retry limit. 9
    s_a = LATIM; %ATIM packet size (slots) 5
    w_a = atim; %ATIM window size (slots) 100

    nn = n; %number of total nodes, including myself 40

    % W_a=4; %contention window size 32
    % m_a=3; %the CW max. 5
    % l_a=3; %the retry limit. 9
    % s_a = 2; %ATIM packet size (slots) 5
    % w_a = 30; %ATIM window size (slots) 100
    % 
    % nn = 5; %number of total nodes, including myself 40


    % w0 = W; %contention window size (slots)
    fprintf('# nodes: %d\n', nn);
    fprintf('CW: %d, max: %d, limit: %d\n', W_a, 2^(log2(W_a)+m_a), l_a-m_a);
    fprintf('====ATIM window related:====\n');
    fprintf('packet size: %d, window size: %d\n', s_a, w_a);


    %for DATA window
    %=========================================
    w_d=data;%100; %DATA window size (slots)
    s_d=LDATA;%30; %DATA packet size (slots)
    W_d=W_a;% con
    m_d=m_a; %m is the CW max.
    l_d=l_a; %j is the retry limit.

    %=========================================

    L_s=36;
    L_c=38;
    L_sc=37;
    Payload = L_s - 16;

    fprintf('====DATA window related:====\n');
    fprintf('packet size: %d, window size: %d\n', s_d, w_d);

    fprintf('lambda: %d\n', lambda);
    fprintf('mu: < %f\n', 1000/((w_a+w_d)/100));


    zw_disp('******************************ATIM start************************************');
    global SYMBOL_VERIFICATION;
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
    % 
    % 
    global atim_p_save_container;
    global atim_pe_container;
    global atim_pf_container;
    global atim_ps_container;
    global atim_sumstate_container;
    global data_p_save_container;    
    global data_pe_container;
    global data_pf_container;
    global data_ps_container;
    global data_sumstate_container;


    global atim_wait_container;
    global data_wait_container;
    global throughput_container_temp;
    global comb_container;
    global atim_energy_container_temp;
    global data_energy_container_temp;

    % zw=1;
    % p=0.5;
    % for i=1:1024
    %     
    %     zw=zw*i*p^i;
    %     
    %     zw_disp(zw);
    % end

    % global p1 p2 p3 p4;


    global WINDOW_ATIM;
    global WINDOW_DATA;


    WINDOW_ATIM = 0;
    WINDOW_DATA = 1;


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



        for i=1:l_a+1
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
    % %        zw_disp(atim_backoff_states_container{1,i});
    % %        zw_disp(atim_windowends_states_container{1,i});
        end

        a=w_d;%20; %DATA window size (slots)
        s=s_d;%4; %DATA packet size (slots)
        W=W_d;
        l=l_d;
        m=m_d;

        for i=1:l_d+1
            if i<m+1
                if SYMBOL_VERIFICATION
                    eval(['C',num2str(i),'=sym(zeros(a*nn,W*2^(i-1)));']);
                    eval(['D',num2str(i),'=sym(zeros(s*nn+1,W*2^(i-1)));']);
                else
                    eval(['C',num2str(i),'=zeros(a*nn,W*2^(i-1));']);
                    eval(['D',num2str(i),'=zeros(s*nn+1,W*2^(i-1));']);
                end
            else
                if SYMBOL_VERIFICATION
                    eval(['C',num2str(i),'=sym(zeros(a*nn,W*2^m));']);
                    eval(['D',num2str(i),'=sym(zeros(s*nn+1,W*2^m));']);
                else
                    eval(['C',num2str(i),'=zeros(a*nn,W*2^m);']);
                    eval(['D',num2str(i),'=zeros(s*nn+1,W*2^m);']);
                end
            end
    % %        zw_disp(atim_backoff_states_container{1,i});
    % %        zw_disp(atim_windowends_states_container{1,i});
        end    

    %     a=w_a;%20; %DATA window size (slots)
    %     s=s_a;%4; %DATA packet size (slots)
    %     W=W_a;
    %     l=l_a;
    %     m=m_a;    

    %     atim_success_probs = zeros(nn,1); 
        %store the successful transmission probabilities for different number of nodes, from the first to the last line, are the prob. of
        %one to nn nodes transmit successfully, respectively.

    %     sum_atim_backoff_states = 0;
    %     sum_atim_windowends_states = 0;
    %     sum_wait_atimwin_ends = 0;

    %     for i=1:nn %generate all A's (i.e. A1 A2 ...) and B's
    %         %A's are for storing the backingoff/freezing states in the ATIM Markov Chain (not including the ATIM window ends state)
    %         %B's are storing the ATIM/DATA window ends state in the ATIM Markov Chain
    %         eval(['global A',num2str(i),';']);
    %         eval(['global B',num2str(i),';']);
    %     end   

        name = sprintf('E:\\save\\atim_mat_%d_%d_%d_%d_%d_%d_%d_%d.mat', nn, W_a, m_a, l_a, s_a, w_a, s_d, w_d);

        try
           load(name);
        catch

            for zw=1:nn%1%
                build_atim_matrix( zw );%this is for obtaining ps, pe, pf of the ATIM stage
            end

            save(name, 'atim_p_save_container', 'atim_pe_container', 'atim_pf_container', ...
                'atim_ps_container', 'atim_sumstate_container',  'atim_energy_container_temp');

        end


        name = sprintf('E:\\save\\atim_wait_%d_%d_%d_%d_%d_%d_%d_%d.mat', nn, W_a, m_a, l_a, s_a, w_a, s_d, w_d);

        try
           load(name);
        catch

            for zw=1:nn%1%
                build_waiting_states( WINDOW_ATIM,  zw );
            end
            save(name, 'atim_wait_container');

        end

        name = sprintf('E:\\save\\data_mat_%d_%d_%d_%d_%d_%d_%d_%d.mat', nn, W_a, m_a, l_a, s_a, w_a, s_d, w_d);

        try
           load(name);
        catch

            for zw=1:nn%1%
                build_data_matrix( zw );
            end
            save(name, 'data_p_save_container', 'data_pe_container', 'data_pf_container', ...
                'data_ps_container', 'data_sumstate_container', 'throughput_container_temp', 'data_energy_container_temp');


        end

        name = sprintf('E:\\save\\data_wait_%d_%d_%d_%d_%d_%d_%d_%d.mat', nn, W_a, m_a, l_a, s_a, w_a, s_d, w_d);

        try
           load(name);
        catch

            for zw=1:nn%1%
                build_waiting_states( WINDOW_DATA,  zw );
            end
            save(name, 'data_wait_container');

        end

        for zw=1:nn%nn%
            throughput_container(zw) = compute_thr( zw );
        end
        
        for zw=1:nn%nn%
            energy_container(zw) = compute_energy( zw );
        end
%         disp(energy_container);
%             return;
        
        for zw=1:nn%nn%
            comb_container{zw} = compute_comb( zw );
        end

    %     return;

        if TEST==1
            t=tic;
            p1 = zeros(nn,nn); 
            p2 = zeros(nn,nn); 
            p3 = zeros(nn,nn);
            p4 = zeros(nn,nn);
            p5 = zeros(nn,nn);
            p6 = zeros(nn,nn);

            all = 0;
            thr(1:nn) = 0;

            for zw=1:nn%
                zw_printf('=======Begin with ATIM for active nodes: %d=======\n', nn-zw+1);

        %         atim_windowends_probs = zeros(nn,1); 
                %store the window ends probability for different number of nodes, from the first to the last line, are the prob. of one to nn nodes
                %not transmit successfully, i.e. the prob. of one to nn nodes left due to window ends.
    %             disp(0.01/nn);
                [pe, pf, ps, sum_temp, ~] = handle_atimwin( zw, 0.01/nn );

                all = all + sum_temp;
                zw_printf('pe:\n');
                zw_disp(pe);
                zw_printf('pf:\n');
                zw_disp(pf);
                DEBUG=0;
                zw_printf('ps:\n');
                zw_disp(ps);
                DEBUG=0;
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
                    [p2e_temp, p2f_temp, p2s_temp, sum_temp, thr_temp, ~] = handle_datawin( zw2, ps(nn+1-zw, nn+1-zw2) );
                    p2e = p2e + p2e_temp;
                    p2f = p2f + p2f_temp;
                    p2s = p2s + p2s_temp;
                    thr(zw) = thr(zw) + thr_temp;
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
    %             DEBUG = 1;
                for zw2=1:nn
                    [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, pef(nn+1-zw, nn+1-zw2) );
    %                 [pew_temp, sum_temp] = waiting_states( WINDOW_DATA, zw2, 1 );
                    pew = pew + pew_temp;
                    all = all + sum_temp;

    %             zw_printf('pef:\n');
    %             zw_disp(pef(nn+1-zw, nn+1-zw2));
    %             zw_printf('pew_temp:\n');
    %             zw_disp(pew_temp);
    %             zw_printf('sum_temp:\n');
    %             zw_disp(sum_temp);
    % 
    %     return;

                end



        % return;

                temp_matrix = zeros(nn,nn);
                temp_matrix(nn+1-zw,:) = sum(pew);
                pew = temp_matrix;

                DEBUG=0;
                zw_printf('pew:\n');
                zw_disp(pew);
                zw_printf('pxn:\n');
                zw_disp(pxn);
                zw_printf('pew*pxn:\n');
                zw_disp(pew*pxn);
                DEBUG=0;
    %

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

    %         disp(thr);

    %         return;
            toc(t);
            DEBUG=0;
            zw_printf('p1:\n');
            zw_disp(p1);
            zw_printf('p2:\n');
            zw_disp(p2);
            zw_printf('p5:\n');
            zw_disp(p5);
            DEBUG=0;
            zw_printf('p3:\n');
            zw_disp(p3);
            zw_printf('p4:\n');
            zw_disp(p4);
    %     return;

        %     ones(1,nn)*p1
        %     ones(1,nn)*p2
        %     ones(1,nn)*p3

            b=sum(p4);
            zw_printf('input for waiting queue...\n');
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
            DEBUG=1;
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
            zw_printf('all: %g, P: %g\n', all, P);



    %         p125=compute_p125(P);
    %         DEBUG=0;
    %         zw_printf('p125:\n');
    %         zw_disp(p125);
    %         L=100;
    %         Ls=(p125-eye(nn))\(-L*ones(nn,1))
    %         
    %         mu=sum(p1+p2+p3+p5+p6)/sum(sum(p1+p2+p3+p5+p6))*p125*Ls+L
    %         lambda=8;
    %         new_P = 1-lambda/(1000/mu)
    %         
    % 
    %         DEBUG=0;


            return;
        end


        time_tic = tic;
        xguess = [];
        for i=1:nn
            xguess = [xguess;0.0005/nn];
        end
        if TEST==0
            xguess = [xguess;P];
        end

    %     for i=1:nn
    %         xguess = [xguess;0.05];
    %     end    

    %     xguess=[1/4 1/4 1/4 1/4 0.01 0.01 0.01 0.01]'; %init and taua

    %     options = optimset('MaxFunEvals',2000, 'MaxIter', 2000);
        diary off;
        options = optimoptions('fsolve','Display','iter');
        xvect = fsolve('solver4d_unsat_math', xguess, options);
        diary on;
        disp(xvect);
    %     thrx=0;
    %     for i=1:nn
    %         thrx = thrx + throughput_container(i)*xvect(i)*nn;%(nn+1-i);
    %     end
        thrx=throughput_container*xvect(1:nn)*nn*54*ceil(PAYLOAD*8/54/10)/LDATA;
        energy=energy_container*xvect(1:nn)*(w_a+w_d)/100;
        fprintf('pernode awake time: %g ms\n', energy);
        
%         if TEST==2
            p125=compute_p125(P);
            DEBUG=0;
            zw_printf('p125:\n');
            zw_disp(p125);
            DEBUG=0;
            L=(w_a+w_d)/100;
            Ls=(p125-eye(nn))\(-L*ones(nn,1));

            mu=((xvect(nn:-1:1))/(sum(xvect(1:nn))))'*p125*Ls+L;            
            lambda = (1-P)*1000/mu;
%         end
%         fprintf('input system throughput: %f Mbps', nn*lambda*PAYLOAD*8/1000000);
        fprintf('actual system throughput: %f Mbps\n', thrx);
%         fprintf('p125:\n');
%         disp(p125);
        L=(w_a+w_d)/100;
%         mu=((xvect(1:nn))/(sum(xvect(1:nn))))'*p125*Ls+L;

        Lt = (Ls-L)*L*2+L*(L-1);

    %     fprintf('Lt:\n');
    %     disp(Lt);

        Ds=(p125-eye(nn))\(-1*Lt);

        zw_printf('Ds+Ls-Ls.^2:\n');
        zw_disp(Ds+Ls-Ls.^2);

        delta_s = ((xvect(nn:-1:1))/(sum(xvect(1:nn))))'.^2*(Ds+Ls-Ls.^2)/10^6;
%         delta_s = ((xvect(nn:-1:1))/(sum(xvect(1:nn))))'*(Ds+Ls-Ls.^2)/10^6;
        fprintf('lambda: %f, mu: %f, service time: %f\n', lambda, 1000/mu, mu);
        
        fprintf('delta_s: %f\n', delta_s);
        mu = 1000/mu;
%         delta_a = 1/lambda;
%         theta = lambda-mu;
%         delta = lambda^3*delta_a^2+lambda*mu^2*delta_s;
% 
%         fprintf('theta: %f\n', theta);
%         fprintf('delta: %f\n', delta);
%         x = [0:1:100];
%         y=1-exp(2*theta/delta*x);
%         plot(y);
% 
%         fun = @(x) exp(2*theta/delta*x);
%         q_len = integral(fun,0,Inf);
%         fprintf('q_len by gg1: %f\n', q_len);

        rho = lambda/mu;
        q_len = rho+(rho^2+lambda^2*delta_s)/(2*(1-rho));
%         q_len = (rho^2+lambda^2*delta_s)/(2*(1-rho));
        fprintf('q_len by mg1: %f\n', q_len);
        fprintf('total delay: %f\n', q_len/lambda);
        toc(time_tic);
        diary off;
        return;
        tau = xvect(1);
        tau2 = xvect(2);
        init = xvect(3);

    %     tau = xvect(3);

        zw_disp('======================================================================');
        zw_disp(['  N = ', num2str(N) ])

        zw_disp('The roots from the default "fsolve" are: ')
    %     zw_disp(['  p = ', num2str(p) ]);
        zw_printf('  init: %g\n', init);

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

        zw_disp(['  tau = ', num2str(tau) ]);
        zw_disp(['  tau2 = ', num2str(tau2) ]);

        taus(runNo)=tau;
    %     Pcs(runNo)=p;

        aggrThr_w(runNo)=N*tau*(1-tau)^(N-1)*Payload/((1-tau)^N+N*tau*(1-tau)^(N-1)*L_s+(1-(1-tau)^N-N*tau*(1-tau)^(N-1))*L_c);
        zw_disp(['  throughput = ', num2str(aggrThr_w(runNo)) ]);


        xaxis(runNo)=N;

        runNo = runNo+1;
    end

    subplot(221);plot(xaxis, taus,'DisplayName','averageTau','YDataSource','averageTau');xlabel('Number of nodes');ylabel('\tau');axis([0,N,0,0.1]);hold all;%figure(gcf);
    subplot(222);plot(xaxis, Pcs,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('Pc');axis([0,N,0,1]);hold all;%figure(gcf);
    subplot(223);plot(xaxis, aggrThr_w,'DisplayName','averageColProb','YDataSource','averageColProb');xlabel('Number of nodes');ylabel('throughput');axis([0,N,0,1]);hold all;%figure(gcf);

end

