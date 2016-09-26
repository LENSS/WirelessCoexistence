function [ thr ] = compute_thr( zw )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global nn;
    [pe, pf, ps, sum_temp] = handle_atimwin( zw, 1 );

    zw_printf('ps:\n');
    zw_disp(ps);
    
    thr = 0;


    for zw2=1:nn
%             [p2e_temp, p2f_temp, p2s_temp, sum_temp] = handle_datawin( zw2, 0.21 );
        [p2e_temp, p2f_temp, p2s_temp, sum_temp, thr_temp] = handle_datawin( zw2, ps(nn+1-zw, nn+1-zw2) );
        thr = thr + thr_temp;

    end


end

