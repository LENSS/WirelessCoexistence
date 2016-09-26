function [ Pe_B, Th_W, ene_B, ene_W ] = compute( sim_time, N_W, N_B, node_struct )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        Pe_B = 1-(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B/2)/(sum([node_struct(N_W/2+1:N_W/2+N_B).success]) + N_B + sum([node_struct(N_W/2+1:N_W/2+N_B).drop]));

        fprintf('Pe_B(%d,%d): %g\n\n', N_B/2, N_W/2, Pe_B);
        fprintf('Th_B(%d,%d): %g\n\n', N_B/2, N_W/2, 1-Pe_B);

        Th_W = (sum([node_struct(1:N_W/2).success]))/(N_W/2)/(sim_time/1e3);
        fprintf('Th_W(%d,%d): %g p/s\n\n', N_B/2, N_W/2, Th_W);

        ene_B = sum([node_struct(N_W/2+1:N_W/2+N_B).awaketime])/(N_B/2)/100/sim_time;
        fprintf('Ene_B(%d,%d): %g\n\n', N_B/2, N_W/2, ene_B);

        ene_W = sum([node_struct(1:N_W/2).awaketime])/(N_W/2)/100/sim_time;
        fprintf('Ene_W(%d,%d): %g\n\n', N_B/2, N_W/2, ene_W);


end

