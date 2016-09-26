function f = servtime_solver( xvect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



    global Sr_B Pe_B Sr_W P_Wa;
    global NN_B NN_W;
    global TR_B TR_W;
%     global T_B T_W;
    global q0b q0w;
    global CASE;

    if CASE==1

        Q0_B=xvect(1);
        Q0_W=xvect(2);

        SS_B = 0;
        SS_W = 0;

        for i=0:NN_B
            beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
        end

        for j=0:NN_W
            beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
        end

    %     for j=0:NN_W+1
    %         SS_B = SS_B + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %         SS_W = SS_W + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %     end

        temp = beta_W(1:NN_W+1) * Sr_B(1:NN_W+1,1:NN_B+1);
        temp1 = beta_B(1:NN_B+1) * Sr_W(1:NN_W+1,1:NN_B+1)';
        temp2 = beta_B(1:NN_B+1) * P_Wa(1:NN_W+1,1:NN_B+1)';

        SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
        SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
        P_W = beta_W(2:NN_W+1) * temp2(2:NN_W+1)';
%         SS_W2 = T_W/P_W;


        S_B = SS_B/(1-beta_B(1))/1e6;
        S_W = SS_W/(1-beta_W(1))/1e6;
%         S_W2 = SS_W2/(1-beta_W(1))/1e6

        q0b = 1-TR_B*S_B;
        q0w = 1-TR_W*S_W;
    %     E_Service = SS/1e6;
    %     traffic_rate_bmac

        f(1) = 1-TR_B*S_B - Q0_B;
        f(2) = 1-TR_W*S_W - Q0_W;
        
    elseif CASE==2
        
        Q0_B=xvect(1);
        
        Q0_W = 0;

        SS_B = 0;
        SS_W = 0;

        for i=0:NN_B
            beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
        end
        
%         disp(beta_B);

        for j=0:NN_W
            beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
        end

    %     for j=0:NN_W+1
    %         SS_B = SS_B + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %         SS_W = SS_W + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %     end

        temp = beta_W(1:NN_W+1) * Sr_B(1:NN_W+1,1:NN_B+1);
        temp1 = beta_B(1:NN_B+1) * Sr_W(1:NN_W+1,1:NN_B+1)';
        temp2 = beta_B(1:NN_B+1) * P_Wa(1:NN_W+1,1:NN_B+1)';

        SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
        SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
        P_W = beta_W(2:NN_W+1) * temp2(2:NN_W+1)';
%         SS_W2 = T_W/P_W;


        S_B = SS_B/(1-beta_B(1))/1e6;
        S_W = SS_W/(1-beta_W(1))/1e6;
%         S_W2 = SS_W2/(1-beta_W(1))/1e6;

        q0b = 1-TR_B*S_B;
        q0w = 1-TR_W*S_W;
    %     E_Service = SS/1e6;
    %     traffic_rate_bmac

        f(1) = 1-TR_B*S_B - Q0_B;
        
    elseif CASE==3
        Q0_W=xvect(1);
        
        Q0_B = 0;

        SS_B = 0;
        SS_W = 0;

        for i=0:NN_B
            beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
        end
        
%         disp(beta_B);

        for j=0:NN_W
            beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
        end

    %     for j=0:NN_W+1
    %         SS_B = SS_B + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %         SS_W = SS_W + nchoosek(N,i)*(1-Q0_B)^i*Q0_B^(N-i)*Service(i);
    %     end

        temp = beta_W(1:NN_W+1) * Sr_B(1:NN_W+1,1:NN_B+1);
        temp1 = beta_B(1:NN_B+1) * Sr_W(1:NN_W+1,1:NN_B+1)';
        temp2 = beta_B(1:NN_B+1) * P_Wa(1:NN_W+1,1:NN_B+1)';

        SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
        SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
        P_W = beta_W(2:NN_W+1) * temp2(2:NN_W+1)';
%         SS_W2 = T_W/P_W;


        S_B = SS_B/(1-beta_B(1))/1e6;
        S_W = SS_W/(1-beta_W(1))/1e6;
%         S_W2 = SS_W2/(1-beta_W(1))/1e6;

        q0b = 1-TR_B*S_B;
        q0w = 1-TR_W*S_W;
    %     E_Service = SS/1e6;
    %     traffic_rate_bmac

        f(1) = 1-TR_W*S_W - Q0_W;
        
    else
        
    end

end

