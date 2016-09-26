function f = main_solver( xvect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

globals;

P_free = xvect(1);
Pe = xvect(2);

A = (Pe+(1-Pe)*(1-tau));


%method 1
% if method == 1
%     Ps = P_free*tau*(1-tau)^(N-1)*P_awake;
% 
%     Pf1 = P_free*tau*(1-tau)^(N-1)*(1-P_awake);
% 
%     Pf2 = P_free*tau*(1-(1-tau)^(N-1));
% 
%     Pns = 1 - P_free*tau;
% 
% else
%method 2
    Ps = P_free*tau*A^(N-1)*P_awake;

    Pf1 = P_free*tau*A^(N-1)*(1-P_awake);

    Pf2 = P_free*tau*(1 - A^(N-1));

    Pns = 1 - P_free*tau;

% end

% P = Ps + Pf1 + Pf2 + Pns;
Pf = 1 - Ps;



Lf = 1/Pf*(Pns*Dns + Pf1*Df1 + Pf2*Df2);
% Lf = Pns*Dns + Pf1*Df1 + Pf2*Df2;

% E_trial = Pf/Ps; %Mean of Geo(Ps)
% V_trial = Pf/Ps^2;%Variance of Geo(Ps)

% T/Lf

% Pe = (Pf)^(T/Lf);

B = Ds/Dns;

f(1) = 1/(1+B*(1-A^N)) - P_free;
f(2) = (Pf)^(T/Lf+1) - Pe;



end

