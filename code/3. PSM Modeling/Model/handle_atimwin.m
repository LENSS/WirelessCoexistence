function [ pe, pf, ps, sum_states, engy ] = handle_atimwin( zw, init )

    global atim_pe_container;
    global atim_pf_container;
    global atim_ps_container;
    global atim_sumstate_container;
    global atim_energy_container_temp;

    
    pe = atim_pe_container{zw}*init;
    pf = atim_pf_container{zw}*init;
    ps = atim_ps_container{zw}*init;
    sum_states = atim_sumstate_container{zw}*init;
    engy = atim_energy_container_temp{zw}*init;

end

