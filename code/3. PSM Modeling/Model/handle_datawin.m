function [ p2e, p2f, p2s, sum_states, thr, engy ] = handle_datawin( zw, init )

    global data_pe_container;
    global data_pf_container;
    global data_ps_container;
    global data_sumstate_container;
    global throughput_container_temp;
    global data_energy_container_temp;

    
    p2e = data_pe_container{zw}*init;
    p2f = data_pf_container{zw}*init;
    p2s = data_ps_container{zw}*init;
    thr = throughput_container_temp{zw}*init;
    engy = data_energy_container_temp{zw}*init;
    sum_states = data_sumstate_container{zw}*init;

end


