function [ p125 ] = compute_p125( P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global DEBUG;
    global nn;
%     global WINDOW_ATIM;
    global WINDOW_DATA;
    
    p1 = zeros(nn,nn); 
    p2 = zeros(nn,nn); 
%     p3 = zeros(nn,nn);
%     p4 = zeros(nn,nn);
    p5 = zeros(nn,nn);

    all = 0;

    for zw=1:nn%
        DEBUG=0;
        zw_printf('=======Begin with ATIM for active nodes: %d=======\n', nn-zw+1);

%         atim_windowends_probs = zeros(nn,1); 
        %store the window ends probability for different number of nodes, from the first to the last line, are the prob. of one to nn nodes
        %not transmit successfully, i.e. the prob. of one to nn nodes left due to window ends.

        [pe, pf, ps, sum_temp] = handle_atimwin( zw, 1 );

        all = all + sum_temp;
        zw_printf('pe:\n');
        zw_disp(pe);
        zw_printf('pf:\n');
        zw_disp(pf);
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
            [p2e_temp, p2f_temp, p2s_temp, sum_temp] = handle_datawin( zw2, ps(nn+1-zw, nn+1-zw2) );
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
%         p3 = p3 + (p2s)*(1-P)*pxn;
%         p4 = p4 + (p2s)*(P)*pxn;

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

    p125=p1+p2+p5;

end

