
clear all;


% syms Pc i W0 j m Lcw Lsw Pf timeslot;
% 
% Lsw + Lcw*(1-Pc)*symsum(Pc^i*i,i, 0, 10) + (1-Pc)/2 * (timeslot+Pf/(1-Pf)) * ( 2*W0*symsum((2*Pc)^i,i,0,m) + ...,
% (2^(m+1))*W0*symsum(Pc^i,i, m+1, 10) - W0*symsum(Pc^i,i,0,10) - symsum(Pc^i*i,i, 0, 10) - symsum(Pc^i,i, 0, 10)  )




%mean expression for 802.11
% syms Pc i W0 j m Lcw Lsw Pf timeslot;
% symsum((2^j*W0-1),j, 0, i) % = W0*(2*2^i - 1) - i - 1
% Lsw + Lcw*(1-Pc)*symsum(Pc^i*i,i, 0, Inf) + (1-Pc)/2 * (timeslot+Pf/(1-Pf)) * ( 2*W0*symsum((2*Pc)^i,i,0,m) + ...,
% (2^(m+1))*W0*symsum(Pc^i,i, m+1, Inf) - W0*symsum(Pc^i,i,0,Inf) - symsum(Pc^i*i,i, 0, Inf) - symsum(Pc^i,i, 0, Inf)  )
% 
% 
% Lsw + Lcw*(1-Pc)*symsum(Pc^i*i,i, 0, Inf) + (1-Pc)/2 * (timeslot+Pf/(1-Pf)) * ( symsum(2*W0*(2*Pc)^i,i,0,m) + ...,
% symsum((2^(m+1))*W0*Pc^i,i, m+1, Inf) - symsum(W0*Pc^i,i,0,Inf) - symsum(Pc^i*i,i, 0, Inf) - symsum(Pc^i,i, 0, Inf)  )


% %variance expression for 802.11
% syms Pc i W0 j m Lcw Pf timeslot;
% 
% syms Var_X Exp_X;
% syms Exp_F Var_F;
% syms Exp_U Var_U;
% syms Pisw Pisb Picw Picb Picbw;
% syms Lsw Lcw Lbt Lmax;
% 
% % Exp_X = Lcw*(1-Pc)*symsum(Pc^i*i,i, 0, Inf) + (1-Pc)/2 * (timeslot+Pf/(1-Pf)) * ( 2*W0*symsum((2*Pc)^i,i,0,m) - ...,
% %   W0*symsum(Pc^i,i,0,Inf) - symsum(Pc^i*i,i, 0, Inf) - symsum(Pc^i,i, 0, Inf) + (2^(m+1))*W0*symsum(Pc^i,i, m+1, Inf) )
% 
% Exp_X = - (Lcw*Pc)/(Pc - 1) - (Pc/2 - 1/2)*(timeslot - Pf/(Pf - 1))*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2^(m + 1)*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)))
%  
% 
% Exp_F = Pf/(1-Pf);
% % Var_F = Pisw*((Lsw-Exp_F)^2)+Pisb*((3*Lbt-Exp_F)^2)+Picw*((Lcw-Exp_F)^2)+Picb*((3*Lbt-Exp_F)^2)+Picbw*((max(Lcw,3*Lbt)-Exp_F)^2);
% Var_F = Pisw*((Lsw-Exp_F)^2)+Pisb*((3*Lbt-Exp_F)^2)+Picw*((Lcw-Exp_F)^2)+Picb*((3*Lbt-Exp_F)^2)+Picbw*((Lmax-Exp_F)^2);
% disp(Var_F);
% 
% Exp_U = (2^j*W0-1)/2;
% Var_U = (4^j*W0^2-1);
% 
% symsum(Exp_U,j,0,i)
% symsum(Var_U,j,0,i)
% ((4*4^i)/3 - 1/3)*W0^2 - i - 1
% 
% % (1-Pc)/2 * Var_F *( 2*W0*symsum((2*Pc)^i,i,0,m) - W0*symsum(Pc^i,i,0,Inf) - symsum(Pc^i*i,i, 0, Inf) - ...,
% %     symsum(Pc^i,i, 0, Inf) + (2*2^m)*W0*symsum(Pc^i,i, m+1, Inf) )
% % 
% % 
% % (1-Pc)/12 * (timeslot+Exp_F)^2 * (4*W0^2/3 * symsum((4*Pc)^i,i,0,m) + (4^(m+1)*W0^2/3)*W0*symsum(Pc^i,i, m+1, Inf) - ...,
% %     W0^2/3*symsum(Pc^i,i,0,Inf) - symsum(Pc^i*i,i, 0, Inf) - symsum(Pc^i,i, 0, Inf))
% 
% 
% -(Pc/2 - 1/2)*(Picw*(Lcw + Pf/(Pf - 1))^2 + Picbw*(Lmax + Pf/(Pf - 1))^2 + Pisw*(Lsw + Pf/(Pf - 1))^2 + ...,
%      Picb*(3*Lbt + Pf/(Pf - 1))^2 + Pisb*(3*Lbt + Pf/(Pf - 1))^2)*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
%      (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2*2^m*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))) ...,
% -(Pc/12 - 1/12)*(timeslot - Pf/(Pf - 1))^2*(1/(Pc - 1) + W0^2/(3*(Pc - 1)) - Pc/(Pc - 1)^2 - (4^(m + 1)*W0^3*(1/(Pc - 1) + ...,
%     (Pc*Pc^m - 1)/(Pc - 1)))/3 + (4*W0^2*(4*4^m*Pc*Pc^m - 1))/(3*(4*Pc - 1)))
%  
% % - (Pc/2 - 1/2)*(Picw*(Lcw + Pf/(Pf - 1))^2 + Picbw*(Lmax + Pf/(Pf - 1))^2 + Pisw*(Lsw + Pf/(Pf - 1))^2 + Picb*(3*Lbt + ...,
% % Pf/(Pf - 1))^2 + Pisb*(3*Lbt + Pf/(Pf - 1))^2)*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - ...,
% % 2*2^m*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))) - (Pc/12 - 1/12)*(timeslot - Pf/(Pf - 1))^2*(1/(Pc - 1) + W0^2/(3*Pc - 3) - ...,
% % Pc/(Pc - 1)^2 - (4^(m + 1)*W0^3*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1)))/3 + (4*W0^2*(4*4^m*Pc*Pc^m - 1))/(12*Pc - 3))
% 
% 
% % (W0*(2*2^i - 1))/2 - i/2 - 1/2
% 
% 
% syms a b c;
% 
% % a = timeslot+Exp_F;
% % b = Lcw;
% % c = Exp_X;
% 
% 
% % (a*((W0*(2*2^i - 1))/2 - i/2 - 1/2)+i*b-c)^2
% 
% % (a^2/4 - a*b + b^2)*i^2 + (a*c - a*b - 2*b*c + (W0*a^2)/2 + a^2/2 - 2^i*W0*a^2 - W0*a*b + 2*2^i*W0*a*b)*i + a*c + (W0*a^2)/2 + a^2/4 + c^2 + (W0^2*a^2)/4 - 2^i*W0*a^2 + W0*a*c - 2^i*W0^2*a^2 + 2^(2*i)*W0^2*a^2 - 2*2^i*W0*a*c
% % 
% % (1-Pc)*symsum(Pc^i*(- 2^i*W0*a^2 + 2*2^i*W0*a*b)*i,i,0,m)
% % 
% % (1-Pc)*symsum(Pc^i*(- 2^m*W0*a^2 + 2*2^m*W0*a*b)*i,i,m+1,Inf)
% % 
% % 
% % (1-Pc)*symsum(Pc^i* (- 2^i*W0*a^2 - 2^i*W0^2*a^2 + 4^i*W0^2*a^2 - 2*2^i*W0*a*c),i,0,m)
% % 
% % (1-Pc)*symsum(Pc^i* (- 2^m*W0*a^2 - 2^m*W0^2*a^2 + 4^m*W0^2*a^2 - 2*2^m*W0*a*c),i,m+1,Inf)
% 
% 
% -(Pc/2 - 1/2)*(Picw*(Lcw + Pf/(Pf - 1))^2 + Picbw*(Lmax + Pf/(Pf - 1))^2 + Pisw*(Lsw + Pf/(Pf - 1))^2 + ...,
%      Picb*(3*Lbt + Pf/(Pf - 1))^2 + Pisb*(3*Lbt + Pf/(Pf - 1))^2)*(1/(Pc - 1) - Pc/(Pc - 1)^2 + W0/(Pc - 1) + ...,
%      (2*W0*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - 2*2^m*W0*(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))) ...,
% -(Pc/12 - 1/12)*(timeslot - Pf/(Pf - 1))^2*(1/(Pc - 1) + W0^2/(3*(Pc - 1)) - Pc/(Pc - 1)^2 - (4^(m + 1)*W0^3*(1/(Pc - 1) + ...,
%     (Pc*Pc^m - 1)/(Pc - 1)))/3 + (4*W0^2*(4*4^m*Pc*Pc^m - 1))/(3*(4*Pc - 1))) + ...,
% a*c + (W0*a^2)/2 + a^2/4 + c^2 + (W0^2*a^2)/4 + ((Pc^2 + Pc)*(a^2/4 - a*b + b^2))/(Pc - 1)^2 + W0*a*c + (Pc*(a*b - a*c + 2*b*c - (W0*a^2)/2 - a^2/2 + W0*a*b))/(Pc - 1) + ...,
% ((2*Pc*W0*a*(a - 2*b))/(2*Pc - 1)^2 - (2^(m + 1)*Pc^(m + 1)*W0*a*(a - 2*b)*(m - 2*Pc*m + 1))/(2*Pc - 1)^2)*(Pc - 1) ..., 
%  -(Pc - 1)*((Pc - Pc*Pc^m - Pc*Pc^m*m + Pc^m*Pc^2*m)/(Pc - 1)^2 - Pc/(Pc - 1)^2)*(2^m*W0*a^2 - 2*2^m*W0*b*a) + ...,
%    (Pc - 1)*((W0^2*a^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) - (W0^2*a^2*(4*4^m*Pc*Pc^m - 1))/(4*Pc - 1) + (W0*a^2*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1) + (2*W0*a*c*(2*2^m*Pc*Pc^m - 1))/(2*Pc - 1)) ...,  
%   -(1/(Pc - 1) + (Pc*Pc^m - 1)/(Pc - 1))*(Pc - 1)*(2^m*W0*a^2 + 2^m*W0^2*a^2 - 4^m*W0^2*a^2 + 2*2^m*W0*a*c)
%  
% 
% % (1-Pc)*(a^2/4 - a*b + b^2)*symsum(Pc^i*i^2,i,0,Inf) + (1-Pc)*(a*c - a*b - 2*b*c + (W0*a^2)/2 + a^2/2 - W0*a*b)*symsum(Pc^i*i,i,0,Inf) + (1-Pc)*(a*c + (W0*a^2)/2 + a^2/4 + c^2 + (W0^2*a^2)/4 + W0*a*c)*symsum(Pc^i,i,0,Inf)
% % 
% % a*c + (W0*a^2)/2 + a^2/4 + c^2 + (W0^2*a^2)/4 + ((Pc^2 + Pc)*(a^2/4 - a*b + b^2))/(Pc - 1)^2 + W0*a*c + (Pc*(a*b - a*c + 2*b*c - (W0*a^2)/2 - a^2/2 + W0*a*b))/(Pc - 1)
% 
% % symsum(Pc^i*i^2,i,0,Inf)
% % -(Pc^2 + Pc)/(Pc - 1)^3
% 
% 
% 


%mean expression for BoXMAC

% syms x timeslot i Wi Wc alpha Lbt
% 3*(1-x)*(Wi-1)/2*timeslot*symsum(x^i,i,0,6)+(alpha/x*1+(1-alpha)*alpha/x*2)*timeslot*(1-x)*symsum(x^i*i,i,0,6)+ ...,
% 3*(1-x)*(Wc-1)/2*timeslot*symsum(x^i*i,i,0,6) + 3*Lbt
% 
% % 3*(1-x)*(Wi-1)/2*timeslot*symsum(x^i,i,0,Inf)+(alpha/x*1+(1-alpha)*alpha/x*2)*timeslot*(1-x)*symsum(x^i*i,i,0,Inf)+ ...,
% % 3*(1-x)*(Wc-1)/2*timeslot*symsum(x^i*i,i,0,Inf) + 3*Lbt


%variance expression for BoXMAC

syms x timeslot i Wi Wc alpha Lbt;
syms Var_Cb Exp_Cb;
syms Var_Xb Exp_Xb;

Exp_Xb = 3*(1-x)*(Wi-1)/2*timeslot*symsum(x^i,i,0,Inf)+(alpha/x*1+(1-alpha)*alpha/2/x*2)*timeslot*(1-x)*symsum(x^i*i,i,0,Inf)+ ...,
  3*(1-x)*(Wc-1)/2*timeslot*symsum(x^i*i,i,0,Inf)

% Exp_Xb = 3*(1-x)*(Wi-1)/2*timeslot*symsum(x^i,i,0,6)+(alpha/x*1+(1-alpha)*alpha/x*2)*timeslot*(1-x)*symsum(x^i*i,i,0,6)+ ...,
%   3*(1-x)*(Wc-1)/2*timeslot*symsum(x^i*i,i,0,6)

syms ab bb cb;


Exp_Cb = alpha/x*1*timeslot+(1-alpha)*alpha/2/x*2*timeslot
Var_Cb = alpha/x*(timeslot-Exp_Cb)^2+(1-alpha)*alpha/2/x*(2*timeslot-Exp_Cb)^2

% ab = 3*timeslot*(Wi-1)/2;
% bb = Exp_Cb+3*timeslot*(Wc-1)/2;
% cb = Exp_Xb;


% symsum(3*(1-x)*(Wi-1)/2*timeslot*x^i,i,0,Inf) + symsum(3*(1-x)*(Wc-1)/2*timeslot*x^i,i,1,Inf)+symsum((alpha/x*1+(1-alpha)*alpha/x*2)*timeslot*(1-x)*x^i*i,i,0,Inf)+3*Lbt

Var_Sb = 9*timeslot^2*(1-x)*(Wi^2-1)/12*symsum(x^i,i,0,Inf) + 9*timeslot^2*(1-x)*(Wc^2-1)/12*symsum(x^i*i,i,0,Inf) + (1-x)*Var_Cb*symsum(x^i*i,i,0,Inf)

% symsum(9*timeslot^2*(1-x)*(Wi^2-1)/12*x^i,i,0,Inf) + symsum(9*timeslot^2*(1-x)*(Wc^2-1)/12*x^i*i,i,0,Inf) + symsum((1-x)*Var_Cb*x^i*i,i,0,Inf)
% 
% (3*timeslot^2*(Wi^2 - 1))/4 - (x*((alpha*(timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x - (alpha*(alpha - 1)*(2*timeslot - (alpha*timeslot)/x + (2*alpha*timeslot*(alpha - 1))/x)^2)/x))/(x - 1) - (3*timeslot^2*x*(Wc^2 - 1))/(4*(x - 1))
% 
% symsum((1-x)*x^i*(ab-cb)^2,i,0,Inf) + symsum((1-x)*x^i*i*(ab-cb)*bb,i,0,Inf) + symsum((1-x)*x^i*i^2*bb^2,i,0,Inf)

(ab - cb)^2 + (bb^2*(x^2 + x))/(x - 1)^2 - (bb*x*(ab - cb))/(x - 1)

