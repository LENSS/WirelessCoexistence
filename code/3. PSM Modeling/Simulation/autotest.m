% y = randn(1,1000);
% h = plot(y(1:200));
% ylim([-5 5]);
% n = 2;
% while n+199 <= numel(y)
%    newy = y(n+(1:199));
%    disp(newy);
%    set(h,'ydata',newy);
%    drawnow;
%    n = n+1;
%    return;
% end
% 
% return;

% T = 100;%s
% T = T*100000;
% n = 30;
% tr = 10;
% PAYLOAD = 1500;
% atim = 100;
% data = 2000 - atim - 1;
% 
% % zw_time=clock();
% % outname=['temp-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
% %     num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
% %     '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];
% % diary(outname);
% % diary on;
% 
% wifi_psm(T, n, PAYLOAD, tr, atim, data);
% % md1(T, n, PAYLOAD, tr, atim, data);
% 
% diary off;
% return;


zw_time(1:6)=0;

T = 100;%s
T = T*100000;
PAYLOAD = 1500;
for atim=200:200:400
    data=2000-atim;
    for n=30:-10:30
        zw_time=clock();
        outname=sprintf('mod(%d,%d,%d)-%d-%d-%d-%d-%d-%g.txt', ...
            n, atim, data, zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
        diary(outname);
        diary on;
        for tr=2.5:2.5:50
            diary on;
            wifi_psm(T, n, PAYLOAD, tr, atim, data);
            diary on;
            fprintf('\n\n');
            diary off;
        end
    end
end