
clear all;
clear all;
clear all;


global showed;
showed = 0;

T_W = 150;% ms
Rho_W = 0.03;
T_B = 50;% ms
Rho_B = 0.15;


fun = @tuning;
x0 = [T_W,Rho_W,T_B,Rho_B];
% x0 = [Rho_W,Rho_B];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [100,0.03,50,0.05];
ub = [200,0.10,100,0.15];
% lb = [0.02,0.02];
% ub = [0.2,0.2];
nonlcon = @nlcon;
x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
% x = ga(fun,4,A,b,Aeq,beq,lb,ub,nonlcon);
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

disp(x);
