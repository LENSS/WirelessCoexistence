
function ret = autotest(test, wifino, bmacno)

%coexist(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b)
% clear all;
zw_time(1:6)=0;
% global TotalTrials;
% TotalTrials = 200;

%varying Wi of BMAC.

switch test
    case{1}%varying Wi
        zw_time=clock();
        outname=['sim4math-Wi-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 150, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 68);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 150, 70, 68);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 68);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 68);
       
    case{2}%varying Wc
        zw_time=clock();
        outname=['sim4math-Wc-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 50, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 30, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 50, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 88);

    case{3}%varying Pb
        zw_time=clock();
        outname=['sim4math-Pb-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 48);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 68);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 108);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 48);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 68);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 88);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 108);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 128);
  case{4}%varying W0
        zw_time=clock();
        outname=['sim4math-W0-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 4, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 6, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 4, 10, 500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 6, 10, 500, 310, 70, 128);
    case{5}%varying Wm
        zw_time=clock();
        outname=['sim4math-Wm-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 6, 8, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 6, 9, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 6, 10, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 4, 8, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 4, 9, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 4, 10, 1000, 310, 70, 128);
    case{6}%varying Pw
        zw_time=clock();
        outname=['sim4math-Pw-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 8, 500, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 8, 1000, 310, 70, 128);
        coexist(10000000, 10, 10, 15, 15, 5, 8, 1500, 310, 70, 128);
    case{7}%varying #device
        zw_time=clock();
        outname=['sim4math-No-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
%         coexist(10000000, 5, 5, 5,  5,  5, 10, 1500, 310, 70, 128);
%         coexist(10000000, 5, 5, 10, 10, 5, 10, 1500, 310, 70, 128);
%         coexist(10000000, 5, 5, 15, 15, 5, 10, 1500, 310, 70, 128);
%         coexist(10000000, 5, 5, 20, 20, 5, 10, 1500, 310, 70, 128);
%         coexist(10000000, 5, 5, 25, 25, 5, 10, 1500, 310, 70, 128);
%         coexist(10000000, 5, 5, 30, 30, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 5, 5, 0, 0, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 10, 10, 0, 0, 5, 10, 1500, 310, 70, 128);
        coexist(10000000, 15, 15, 0, 0, 5, 10, 1500, 310, 70, 128);
    case{8}%varying #device
        zw_time=clock();
        outname=['sim4math-No-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
%         coexist(15000000, 10, 10, 5,  5,  5, 10, 1500, 310, 70, 128);
%         coexist(15000000, 10, 10, 10, 10, 5, 10, 1500, 310, 70, 128);
%         coexist(15000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
%         coexist(15000000, 10, 10, 20, 20, 5, 10, 1500, 310, 70, 128);
%         coexist(15000000, 10, 10, 25, 25, 5, 10, 1500, 310, 70, 128);
%         coexist(15000000, 10, 10, 30, 30, 5, 10, 1500, 310, 70, 128);
         coexist(10000000, 0, 0, 10, 10, 5, 10, 1500, 310, 70, 128);
         coexist(10000000, 0, 0, 30, 30, 5, 10, 1500, 310, 70, 128);
    case{9}%varying #device
        zw_time=clock();
        outname=['sim4math-No-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
%         coexist(20000000, 15, 15, 5,  5,  5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 10, 10, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 15, 15, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 20, 20, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 25, 25, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 30, 30, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 5, 5, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 1, 1, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 6, 6, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 0, 0, 10, 10, 5, 10, 1500, 310, 70, 128);
        
        
    otherwise
%          coexist(100000, 0, 0, 30, 60, 0, 0, 0, 160, 60, 100);
        zw_time=clock();
        outname=['sim4math-', num2str(zw_time(2)),'-', num2str(zw_time(3)),'-', ...
            num2str(ceil(zw_time(1))),'-',num2str(zw_time(4)), ...
            '-', num2str(zw_time(5)),'-', num2str(zw_time(6)), '.txt'];

        diary(outname);
%          coexist(1000000, 40, 40, 20, 20, 5, 10, 1500, 160, 60, 100);

%          coexist(1000000, 10, 10, 10, 10, 5, 10, 1500, 310, 70, 128);
%          coexist(1000000, 10, 10, 10, 10, 5, 10, 1500, 160, 60, 128);
%          return;
         
%         coexist(20000000, 15, 15, 5,  5,  5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 10, 10, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 15, 15, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 20, 20, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 25, 25, 5, 10, 1500, 310, 70, 128);
%         coexist(20000000, 15, 15, 30, 30, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 5, 5, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 1, 1, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 6, 6, 0, 0, 5, 10, 1500, 310, 70, 128);
%         coexist(1000000, 0, 0, 10, 10, 5, 10, 1500, 310, 70, 128);
    for Packet_w=500:500:500%[500,1000,1500]
    for aMinBE=4:6
    for aMaxBE=8:10;
    for Packet_b=128:-20:48%[128,108,88,68,48]
    for Wc_b=70:-20:30%[70,50,30]
    for Wi_b=310:-80:70%[310,230,150,70]
        coexist(2000000, wifino, wifino, bmacno, bmacno, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
    end
    end
    end
    end
    end
    end

end


% WIFI_START=0;
% WIFI_END=0;
% 
% BMAC_START=1;
% BMAC_END=1;

% for Packet_w=1500:500:1500%[500,1000,1500]
% for aMinBE=6:6
% for aMaxBE=10:10;
% 
% for Packet_b=128:-20:128%[128,108,88,68,48]
% for Wc_b=70:-20:70%[70,50,30]
% for Wi_b=310:-80:310%[310,230,150,70]

% coexist(500000, 2, 2, 4, 4, 4, 10, 1500, 310, 70, 128);
% coexist(500000, 2, 2, 4, 4, 5, 10, 1500, 310, 70, 128);
% coexist(500000, 2, 2, 4, 4, 6, 10, 1500, 310, 70, 128);
% coexist(500000, 2, 2, 4, 4, 5, 10, 1500, 310, 70, 128);
% coexist(500000, 2, 2, 4, 4, 6, 10, 1500, 310, 70, 128);
ret = 0;

