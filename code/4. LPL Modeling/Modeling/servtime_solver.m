function f = servtime_solver( xvect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    globals;


    global Service;

    P0=xvect(1);
    
    SS = 0;
    for i=1:N
        SS = SS + nchoosek(N,i)*(1-P0)^i*P0^(N-i)*Service(i);
    end
    
    E_Service = SS/(1-nchoosek(N,0)*(1-P0)^0*P0^(N))/1e6;
%     E_Service = SS/1e6;
%     traffic_rate_bmac
    
    f(1) = 1-traffic_rate_bmac*E_Service - P0;

end

