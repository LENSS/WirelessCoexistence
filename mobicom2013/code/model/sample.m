clear all
% disp('Solve the following set of nonlinear algebraic equations:')
% disp('  sin(x) + y^2 + ln(z) - 7 = 0 ')
% disp('  3*x + 2^y - z^3 + 1      = 0 ')
% disp('  x + y + z - 5            = 0 ')
% disp(' ')
m=4;
n=3;
xguess=[0 0 0]';
xvect = fsolve('verify', xguess);
alpha = xvect(1);
beta = xvect(2);
Pcoll = xvect(3);
disp('The roots from the default "fsolve" are: ')
disp(['  alpha = ', num2str(alpha) ])
disp(['  beta = ', num2str(beta) ])
disp(['  Pcoll = ', num2str(Pcoll) ])
x = alpha+(1-alpha)*beta;
y = Pcoll*(1-x^(m+1));
disp(['  x = ', num2str(x) ])
disp(['  y = ', num2str(y) ])
Pcf = (x^(m+1))*(1-y^(n+1))/(1-y);
Pcr = y^(n+1);
R = 1-Pcf-Pcr;
disp(['  Pcf = ', num2str(Pcf) ])
disp(['  Pcr = ', num2str(Pcr) ])
disp(['  R = ', num2str(R) ])
% % Repeat with a different set of options -------------------------------
% options(2) = 1.e-6;        %Tolerance for x
% options(3) = 1.e-6;        %Tolerance for f
% options(5) = 1;            %Levenberg-Marquardt Method
% xvect = fsolve('myfun', xguess, options);
% x = xvect(1);
% y = xvect(2);
% z = xvect(3);
% disp('The roots from "fsolve" with Levenberg-Marquardt are: ')
% disp(['  x = ', num2str(x) ])
% disp(['  y = ', num2str(y) ])
% disp(['  z = ', num2str(z) ])


% H=0.32;
% Pc0=0.23;W=0.18;
% x0 = [2*W; Pc0+2*H]; % ???
% options = optimset('Display','off'); 
% k=0:0.01:1; % ??????[0 1]
% for i=1:1:length(k)
% kk=k(i);
% x = fsolve(@(x) myfun(x,kk), x0, options);%????????
% x1(i)=x(1);
% x2(i)=x(2);
% end
% plot(k,x1,'-b',k,x2,'-r');
% xlabel('k')
% legend('x1','x2')

