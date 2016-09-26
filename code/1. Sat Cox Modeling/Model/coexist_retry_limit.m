% m=4;
% n=3;
xguess=[0 0 0]';
xvect = fsolve('coexist_math_retry_limit', xguess);
alpha = xvect(1);
Pc = xvect(2);
Pf = xvect(3);

disp('The roots from the default "fsolve" are: ')
disp(['  alpha = ', num2str(alpha) ])
disp(['  Pc = ', num2str(Pc) ])
disp(['  Pf = ', num2str(Pf) ])