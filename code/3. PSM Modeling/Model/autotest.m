% y = evrnd(0,3,100,1);
% cdfplot(y)
% 
% return;
clearvars;
clear all;
clear all;
clear all;

% x = [0:0.01:1];
% y=exp(-8*x);
% fun = @(x) exp(-8*x);
% q = integral(fun,0,Inf);
% 
% plot(y);
% 
% return;
global TEST;
TEST = 2;
% 0, no test; 
% 1, single run test
% 2, test with fixed P

%solver4d_unsat( n, payload, tr, atim, data, p0 )
% solver4d_unsat(30, 1500, 2.5, 100, 1900, 0.95);
% solver4d_unsat(30, 1500, 5, 100, 1900, 0.90);
% solver4d_unsat(30, 1500, 7.5, 100, 1900, 0.85);
% solver4d_unsat(30, 1500, 10, 100, 1900, 0.76);
% return;


zw_time(1:6)=0;

T = 10000000;
PAYLOAD = 1500;
for atim=100:100:200%300
    data=2000-atim;
    for n=30:-10:10
        zw_time=clock();
%         outname=['mod-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
%             num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
%             '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];
        outname=sprintf('mod(%d,%d,%d)-%d-%d-%d-%d-%d-%g.txt', ...
            n, atim, data, zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
        diary(outname);
        diary on;
        for tr=5:5:5
            for p0=0.9:-0.1:0
                diary on;
                solver4d_unsat(n, PAYLOAD, tr, atim, data, p0);
                diary on;
                fprintf('\n\n');
                diary off;
            end
        end
    end
end