function [ SS_B, SS_W, Q0_B, en_W ] = compute(NN_B, TR_B, T_B, ...
                                    NN_W, TR_W, T_W)

    global CASE;
    global q0b q0w;
    global Sr_B Pe_B Sr_W En_W P_Wa P_Wd;
    
    
    iter = 'final';%'off';


    if TR_B<(1e6/T_B) && TR_W<(1e6/T_W)
        CASE = 1;

        Q0_B = 0.5;
        Q0_W = 0.5;

        xguess=[Q0_B Q0_W]';
        options = optimoptions('fsolve', 'MaxFunEvals',...
            1000,'TolFun',1e-12,'TolX',1e-12,...
            'Display',iter); %'off' 'final'

        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);

        exitflag;
%         output
        if exitflag<=0

%             ii = ii + 1;
%             result(ii,1) = 0;
%             result(ii,2) = 0;
%             result(ii,3) = 0;
%             result(ii,4) = 0;
%             zzww{ii,1} = sprintf('N_W:%d,TR_W:%d,N_B:%d,TR_B:%g', NN_W, TR_W, NN_B, TR_B);
%            continue;
            Q0_B = q0b;
            Q0_W = q0w;

        else
            q0b;
            q0w;
            Q0_B = xvect(1);
            Q0_W = xvect(2);
        end
%             fprintf('Expected Q0_B: %g\n', Q0_B);
%             fprintf('Expected Q0_W: %g\n', Q0_W);


    elseif TR_B<(1e6/T_B) && TR_W>=(1e6/T_W)
        CASE = 2;
        Q0_B = 0.5;
        Q0_W = 0;

        xguess=[Q0_B]';
        options = optimoptions('fsolve', 'MaxFunEvals',...
            1000,'TolFun',1e-12,'TolX',1e-12,...
            'Display',iter); %'off' 'final'

        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);

        exitflag;

        if exitflag<=0

            Q0_B = q0b;
%                 Q0_W = q0w


        else
            q0b;
%                 q0w
            Q0_B = xvect(1);
        end            
%             fprintf('Expected Q0_B: %g\n', Q0_B);
%             fprintf('Expected Q0_W: %g\n', Q0_W);

    elseif TR_B>=(1e6/T_B) && TR_W<(1e6/T_W)
        CASE = 3;

        Q0_B = 0;
        Q0_W = 0.001;

        xguess=[Q0_W]';
        options = optimoptions('fsolve', 'MaxFunEvals',...
            1000,'TolFun',1e-12,'TolX',1e-12,...
            'Display',iter); %'off' 'final'

        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);

        exitflag;

        if exitflag<=0

            Q0_W = q0w;


        else
%                 q0b
            q0w;
            Q0_W = xvect(1);
        end            
%             fprintf('Expected Q0_B: %g\n', Q0_B);
%             fprintf('Expected Q0_W: %g\n', Q0_W);

    else %if TR_B>=(1e6/T_B) && TR_W>=(1e6/T_W)
        CASE = 4;
        Q0_B = 0;
        Q0_W = 0;

    end

    if Q0_B<0 && Q0_W>0 && Q0_W<1

        CASE = 3;

        Q0_B = 0;
        Q0_W = 0.5;

        xguess=[Q0_W]';
        options = optimoptions('fsolve', 'MaxFunEvals',...
            1000,'TolFun',1e-12,'TolX',1e-12,...
            'Display',iter); %'off' 'final'

        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);

        exitflag;

        if exitflag<=0

%             ii = ii + 1;
%             result(ii,1) = 0;
%             result(ii,2) = 0;
%             result(ii,3) = 0;
%             result(ii,4) = 0;
%             zzww{ii,1} = sprintf('N_W:%d,TR_W:%d,N_B:%d,TR_B:%g', NN_W, TR_W, NN_B, TR_B);
%            continue;
%                 Q0_B = q0b
            Q0_W = q0w;


        else
%                 q0b
            q0w;
            Q0_W = xvect(1);
        end            

%                 error('error2');
% 
%             ii = ii + 1;
% 
%             result_ThrB2{ii,jj} = -2;
%             result_ThrW2{ii,jj} = -2;
%             result_EneB2{ii,jj} = -2;
%             result_EneW2{ii,jj} = -2;
% 
%             result_ThrB2{ii,1} = sprintf('TR_W=%g', TR_W);
%             result_ThrW2{ii,1} = sprintf('TR_W=%g', TR_W);
%             result_EneB2{ii,1} = sprintf('TR_W=%g', TR_W);
%             result_EneW2{ii,1} = sprintf('TR_W=%g', TR_W);
%             continue;

    elseif Q0_B>0 && Q0_W<0 && Q0_B<1

        CASE = 2;
        Q0_B = 0.5;
        Q0_W = 0;

        xguess=[Q0_B]';
        options = optimoptions('fsolve', 'MaxFunEvals',...
            1000,'TolFun',1e-12,'TolX',1e-12,...
            'Display',iter); %'off' 'final'

        [xvect,fval,exitflag,output]  = fsolve('servtime_solver', xguess, options);

        exitflag;

        if exitflag<=0

            Q0_B = q0b;
%                 Q0_W = q0w


        else
            q0b;
%                 q0w
            Q0_B = xvect(1);
        end            

    elseif Q0_B<0 && Q0_W<0
        Q0_B = 0;
        Q0_W = 0;

    elseif Q0_B==0 && (Q0_W<0 || Q0_W>=1)
        Q0_B = 0;
        Q0_W = 0;


    elseif Q0_W==0 && (Q0_B<0 || Q0_B>=1 )
        Q0_B = 0;
        Q0_W = 0;

    end

    fprintf('Expected Q0_B: %g\n', Q0_B);
    fprintf('Expected Q0_W: %g\n', Q0_W);


    for i=0:NN_B
        beta_B(i+1) = nchoosek(NN_B,i)*(1-Q0_B)^i*Q0_B^(NN_B-i);
    end

    for j=0:NN_W
        beta_W(j+1) = nchoosek(NN_W,j)*(1-Q0_W)^j*Q0_W^(NN_W-j);
    end

    beta_W;
    beta_B;
%     Pe_B

    temp = beta_W(1:NN_W+1) * Pe_B(1:NN_W+1,1:NN_B+1);
    temp1 = beta_B(1:NN_B+1) * (1e6./Sr_W(1:NN_W+1,1:NN_B+1))';
    temp2 = beta_B(1:NN_B+1) * (En_W(1:NN_W+1,1:NN_B+1))';
    temp3 = beta_B(1:NN_B+1) * P_Wa(1:NN_W+1,1:NN_B+1)';
    
    (1e6./Sr_W)';
    temp1;

    SS_B = beta_B(2:NN_B+1) * temp(2:NN_B+1)';
    SS_W = beta_W(2:NN_W+1) * temp1(2:NN_W+1)';
    en_W = beta_W(2:NN_W+1) * temp2(2:NN_W+1)';
    P_W = beta_W(2:NN_W+1) * temp3(2:NN_W+1)';
    SS_W2 = 1e6/(T_W/P_W);
    
    SS_W;

    SS_B = SS_B/(1-beta_B(1));
    SS_W = SS_W/(1-beta_W(1));
    if SS_W > TR_W
       SS_W = TR_W; 
    end



    en_W = en_W/(1-beta_W(1));



end

