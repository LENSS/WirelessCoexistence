function [ energy ] = compute_energy( zw )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global nn;
    [pe, pf, ps, sum_temp, atim_engy] = handle_atimwin( zw, 1 );

%     fprintf('ps:\n');
%     disp(ps);
%     fprintf('pe:\n');
%     disp(pe);
%     fprintf('pf:\n');
%     disp(pf);
    
    energy = atim_engy;
%     fprintf('energy: %g\n', energy);
    

    for zw2=1:nn
%             [p2e_temp, p2f_temp, p2s_temp, sum_temp] = handle_datawin( zw2, 0.21 );
        [p2e_temp, p2f_temp, p2s_temp, sum_temp, thr_temp, data_engy] = handle_datawin( zw2, ps(nn+1-zw, nn+1-zw2) );
        energy = energy + data_engy;
%         fprintf('ps: %g, energy: %g\n', ps(nn+1-zw, nn+1-zw2), energy);

    end

%         fprintf('\n');

end

