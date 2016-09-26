function [ pw, sum_states ] = waiting_states( window_type, zw, input )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%     global nn;
%     global w_d;    
%     global w_a;    
%     global s_d;    
%     global s_a;    
%     global data_p_save_container;
%     global atim_p_save_container;
%     global WINDOW_ATIM;
    global WINDOW_DATA;
    global atim_wait_container;
    global data_wait_container;
    
   
    if window_type == WINDOW_DATA
        pw = data_wait_container{1,zw}*input;
        sum_states = data_wait_container{2,zw}*input;
    else
        pw = atim_wait_container{1,zw}*input;
        sum_states = atim_wait_container{2,zw}*input;
    end

end

