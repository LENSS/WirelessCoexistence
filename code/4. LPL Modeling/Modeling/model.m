
function out = model( sim, nodes, timeout, cycle, ratio, cw, traffic  )



% clear;
% clear;
% clear;
% format long;


globals;


tau = 2/3/(cw+10+1);%1/40;

P_awake = ratio;%1/10;

N = nodes;

Ds = 1020;%3000;%us
Df1 = Ds + timeout;%us
Df2 = Ds + timeout;%us
Dns = 10;%us
T = cycle; %1000000;

% tau = 1/40;
% 
% P_awake = 1/10;
% 
% N = 40;
% 
% Ds = 3000;%us
% Df1 = 300;%us
% Df2 = 300;%us
% Dns = 30;%us
% T = 200000;


P_free = 0.01;
Pe = 0.01;

xguess=[P_free Pe]';
options = optimoptions('fsolve','Display','iter'); %'off'
xvect = fsolve('main_solver', xguess, options);
% disp(xvect);

P_free = xvect(1);
Pe = xvect(2);

fprintf('P_free(%d): %g\n', N, P_free);
fprintf('Pe(%d): %g\n', N, Pe);

out(1) = P_free;
out(2) = Pe;

% following is the new idea
%============================================================================================
% 
E_service = (1-Pe)*T/2 + Pe*T;
out(4) = E_service;
fprintf('E_service(%d): %g\n', N, E_service);
% 
%============================================================================================

return;

% following is the newest idea, still not good (due to the uniform distribution of Pr(awake))
%============================================================================================

A = (Pe+(1-Pe)*(1-tau));


%method 1
if method == 1
    Ps = P_free*tau*(1-tau)^(N-1)*P_awake;

    Pf1 = P_free*tau*(1-tau)^(N-1)*(1-P_awake);

    Pf2 = P_free*tau*(1-(1-tau)^(N-1));

    Pns = 1 - P_free*tau;

else
%method 2
    Ps = P_free*tau*A^(N-1)*P_awake;

    Pf1 = P_free*tau*A^(N-1)*(1-P_awake);

    Pf2 = P_free*tau*(1 - A^(N-1));

    Pns = 1 - P_free*tau;

end




% P = Ps + Pf1 + Pf2 + Pns;
Pf = 1 - Ps;


% return;
% B = Ds/Dns;

% S = P_free * N*(1-Pe)*P_awake*tau*A^(N-1)*B;
% out(3) = S;
% fprintf('S(%d): %g\n', N, S);

Lf = 1/Pf*(Pns*Dns + Pf1*Df1 + Pf2*Df2);
E_service = 0;
for i=0:ceil(T/Lf)
%     Ps = Ps*1.0002
%     Pf = 1 - Ps
    E_service = E_service + Pf^i*Ps*(i+1)*Lf;
    
end

E_service = E_service + Pe*T;

% E_trial = Pf/Ps; %Mean of Geo(Ps)
% % V_trial = Pf/Ps^2; %Variance of Geo(Ps)
% 
% E_service = min(E_trial * Lf,T);
% P0 = 1- traffic * E_service/1e6;
out(4) = E_service;
fprintf('E_service(%d): %g\n', N, E_service);

%============================================================================================





% % following is the original idea, which is not correct.
% %============================================================================================
% 
% A = (Pe+(1-Pe)*(1-tau));
% 
% 
% %method 1
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
% %method 2
%     Ps = P_free*tau*A^(N-1)*P_awake;
% 
%     Pf1 = P_free*tau*A^(N-1)*(1-P_awake);
% 
%     Pf2 = P_free*tau*(1 - A^(N-1));
% 
%     Pns = 1 - P_free*tau;
% 
% end
% 
% % P = Ps + Pf1 + Pf2 + Pns;
% Pf = 1 - Ps;
% 
% B = Ds/Dns;
% 
% % S = P_free * N*(1-Pe)*P_awake*tau*A^(N-1)*B;
% % out(3) = S;
% % fprintf('S(%d): %g\n', N, S);
% 
% Lf = 1/Pf*(Pns*Dns + Pf1*Df1 + Pf2*Df2);
% 
% E_trial = Pf/Ps; %Mean of Geo(Ps)
% % V_trial = Pf/Ps^2; %Variance of Geo(Ps)
% 
% E_service = min(E_trial * Lf,T);
% % P0 = 1- traffic * E_service/1e6;
% out(4) = E_service;
% fprintf('E_service(%d): %g\n', N, E_service);
% % 
% % E_service2 = 0;
% % temp = ceil((T-Ds)/Lf);
% % for i=0:temp
% %   E_service2 = E_service2 + (1-Ps)^i*Ps*(i*Lf+Ds);
% % end
% % E_service2 = E_service2 + (1-Ps)^((T-Ds)/Lf+1)*T;
% % 
% % 
% % out(3) = E_service2;
% % % fprintf('E_service2(%d): %g\n', N, E_service2);
% %============================================================================================


% T/Lf

end
