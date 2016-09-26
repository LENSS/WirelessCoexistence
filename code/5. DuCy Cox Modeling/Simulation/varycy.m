function [ output_args ] = varycy(sim_time, N_B, L_B, Ackto_B, T_B, Rho_B, CW_B, TR_B, ...
                            N_W, L_W, Ackto_W, T_W, Rho_W, CW_W, TR_W )

    global packet_mactime;
    global packet_sucinatim;

    T_W=[50 50 100 100 200 200];
    Rho_W=[0.1 0.2 0.05 0.1 0.025 0.05];

    zw_printf('LPL cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
    zw_printf('PSM cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);

    outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...% 500, ...%
                            N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                            N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
                        
%     disp(outname);

    node_struct = cox(sim_time, N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
                                N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);

    save(outname, 'node_struct', 'packet_mactime', 'packet_sucinatim');


%     T_W=400;
%     Rho_W=0.04;
% 
%     zw_printf('LPL cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
%     zw_printf('PSM cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
% 
%     outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...% 500, ...%
%                             N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                             N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
%                         
%     disp(outname);
% 
%     node_struct = cox(sim_time, N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                                 N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
% 
%     save(outname, 'node_struct', 'packet_mactime', 'packet_sucinatim');
% 
% 
%     T_W=800;
%     Rho_W=0.0125;
% 
%     zw_printf('LPL cycle: %d ms, duty-cycle ratio: %g\n', T_B, Rho_B);
%     zw_printf('PSM cycle: %d ms, duty-cycle ratio: %g\n', T_W, Rho_W);
% 
%     outname=sprintf('mat\\stat(%d)(%d-%d-%d-%d-%g-%d-%g)(%d-%d-%d-%d-%g-%d-%g).mat', sim_time/1000, ...% 500, ...%
%                             N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                             N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
%                         
%     disp(outname);
% 
%     node_struct = cox(sim_time, N_B, L_B, Ackto_B, T_B, Rho_B, CW_B(1), TR_B, ...
%                                 N_W, L_W, Ackto_W, T_W, Rho_W, CW_W(1), TR_W);
% 
%     save(outname, 'node_struct', 'packet_mactime', 'packet_sucinatim');



end

