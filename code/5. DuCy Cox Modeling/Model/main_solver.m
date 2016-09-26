function f = main_solver( xvect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

globals;

% Pi = xvect(1);
% Pe_B = xvect(2);
% 
% tau_B = 0;
% Pa_B = 0;
% N_B = 0;
% 
% A = (Pe_B+(1-Pe_B)*(1-tau_B));
% 
% Ps = Pi*tau_B*A^(N_B-1)*Pa_B;
% 
% Pf1 = Pi*tau_B*A^(N_B-1)*(1-Pa_B);
% 
% Pf2 = Pi*tau_B*(1 - A^(N_B-1));
% 
% Pns = 1 - Pi*tau_B;
% 
% Pf = 1 - Ps;
% Lf = 1/Pf*(Pns*Dns + Pf1*Df1 + Pf2*Df2);
% B = Ds/Dns;
% 
% f(1) = 1/(1+B*(1-A^N_B)) - Pi;
% f(2) = (Pf)^(T/Lf+1) - Pe_B;



Pi = xvect(1);
Pe_B = xvect(2);
tau_W = xvect(3);

temp_B = Pe_B+(1-Pe_B)*(1-tau_B);
temp_W = Pe_W+(1-Pe_W)*(1-tau_W)^X_W;



%Channel model
Pib0 = N_B*tau_B*(1-Pe_B)*temp_B^(N_B-1)*Pa_B*temp_W;
Pib1 = N_B*tau_B*(1-Pe_B)*temp_B^(N_B-1)*(1-Pa_B)*temp_W;
Pib2 = (1-temp_B^N_B-N_B*tau_B*(1-Pe_B)*temp_B^(N_B-1))*temp_W;
Pib3 = (1-Pe_W)*X_W*tau_W*(1-tau_W)^(X_W-1)*temp_B^N_B;
Pib4 = ((1-Pe_W)*(1-(1-tau_W)^(X_W)-X_W*tau_W*(1-tau_W)^(X_W-1)))*temp_B^N_B;
Pib5 = (1-temp_B^N_B)*(1-Pe_W)*(1-(1-tau_W)^X_W);


if N_B ~= 0
%LPL model
    Ps_B = Pi*tau_B*temp_B^(N_B-1)*temp_W;
    Pf1_B = 0;
%     Ps_B = Pi*tau_B*temp_B^(N_B-1)*Pa_B*temp_W;
%     Pf1_B = Pi*tau_B*temp_B^(N_B-1)*(1-Pa_B)*temp_W;
    Pf2_B = Pi*tau_B*(1 - temp_B^(N_B-1))*temp_W;
    Pf3_B = Pi*tau_B*(1-Pe_W)*(1-(1-tau_W)^X_W);
    Pi_B = 1 - Pi*tau_B;


    % P = Ps + Pf1 + Pf2 + Pns;
    Pf_B = 1 - Ps_B;


    ALf_B = (Pi_B*Di_B + (Pf1_B + Pf2_B + Pf3_B)*Df_B)/Pf_B;%
    % Lf = Pns*Dns + Pf1*Df1 + Pf2*Df2;

    % E_trial = Pf/Ps; %Mean of Geo(Ps)
    % V_trial = Pf/Ps^2;%Variance of Geo(Ps)

    % T/Lf

    % Pe = (Pf)^(T/Lf);

    Ls_B = Ds_B/timeslot;
    Lf_B = Df_B/timeslot;
    Ls_W = Ds_W/timeslot;
    Lf_W = Df_W/timeslot;
    Lf_BW = Lf_B;


    f(1) = 1/(1+Pib0*Ls_B+Pib1*Lf_B+Pib2*Lf_B+Pib3*Ls_W+Pib4*Lf_W+Pib5*Lf_BW) - Pi;
    f(2) = (Pf_B)^((Ta_B-Df_B/2)/ALf_B+1) - Pe_B;
%     f(2) = (Pf_B)^(T_B/ALf_B+1) - Pe_B;
else
    f(1) = 1 - Pi;
    f(2) = 0 - Pe_B;
end



p = 1-((1-tau_W)^(X_W-1)*temp_B^N_B);
% p = Pib0 + Pib1 + Pib2 + Pib3 + Pib4 + Pib5;


% if X_W > 1
    f(3) = 2*(1-2*p)/((1-2*p)*(W0_W+1)+p*W0_W*(1-(2*p)^m_W)) - tau_W; %PSM model, Bianchi
% else
%     f(3) = (1-2*p)/((1-2*p)*(W0_W+1)+p*W0_W*(1-(2*p)^m_W))/(5-X_W) - tau_W; %PSM model, Bianchi
% end


end

