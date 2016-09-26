
function ret = autotest2(test, wifino, bmacno)

%coexist_channel(T, WIFI_START, WIFI_END, BMAC_START, BMAC_END, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b)
zw_time(1:6)=0;
% digits(32);
warning('off','all');
warning;
%varying Wi of BMAC.

switch test
    case{1}%varying Wi
        zw_time=clock();
        outname='logs/math_Wi2.txt';

        diary(outname);
        diary on;
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 150, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 68);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 150, 70, 68);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 68);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 68);
        diary off;
       
    case{2}%varying Wc
        zw_time=clock();
        outname='logs/math_Wc2.txt';

        diary(outname);
        diary on;
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 50, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 30, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 50, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 70, 70, 88);
        diary off;

    case{3}%varying Pb
        zw_time=clock();
        outname='logs/math_Pb2.txt';

        diary(outname);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 48);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 68);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 108);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 48);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 68);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 88);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 108);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 230, 30, 128);
  case{4}%varying W0
        zw_time=clock();
        outname='logs/math_W02.txt';

        diary(outname);
        coexist_channel(10000000, 10, 10, 15, 15, 4, 10, 1500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 6, 10, 1500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 4, 10, 500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 6, 10, 500, 310, 70, 128);
    case{5}%varying Wm
        zw_time=clock();
        outname='logs/math_Wm2.txt';

        diary(outname);
        coexist_channel(10000000, 10, 10, 15, 15, 6, 8, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 6, 9, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 6, 10, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 4, 8, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 4, 9, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 4, 10, 1000, 310, 70, 128);
    case{6}%varying Pw
        zw_time=clock();
        outname='logs/math_Pw2.txt';

        diary(outname);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 8, 500, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 8, 1000, 310, 70, 128);
        coexist_channel(10000000, 10, 10, 15, 15, 5, 8, 1500, 310, 70, 128);
    case{7}%varying #device
        zw_time=clock();
        outname='logs/math_No1.txt';

        diary(outname);
        coexist_channel(10000000, 5, 30, 0, 0, 5, 10, 1500, 320, 80, 128);
        coexist_channel(10000000, 0, 0, 5, 30, 5, 10, 1500, 320, 80, 128);
        coexist_channel(10000000, 5, 5, 10, 10, 5, 10, 1500, 320, 80, 128);
        coexist_channel(15000000, 10, 10, 10, 10, 5, 10, 1500, 320, 80, 128);
        coexist_channel(20000000, 15, 15, 10, 10, 5, 10, 1500, 320, 80, 128);
        coexist_channel(10000000, 5, 5, 30, 30, 5, 10, 1500, 320, 80, 128);
        coexist_channel(15000000, 10, 10, 30, 30, 5, 10, 1500, 320, 80, 128);
        coexist_channel(20000000, 15, 15, 30, 30, 5, 10, 1500, 320, 80, 128);
    case{8}%varying #device
        zw_time=clock();
        outname='logs/math_No2.txt';

        diary(outname);
        coexist_channel(15000000, 10, 10, 5,  5,  5, 10, 1500, 310, 70, 128);
        coexist_channel(15000000, 10, 10, 10, 10, 5, 10, 1500, 310, 70, 128);
        coexist_channel(15000000, 10, 10, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist_channel(15000000, 10, 10, 20, 20, 5, 10, 1500, 310, 70, 128);
        coexist_channel(15000000, 10, 10, 25, 25, 5, 10, 1500, 310, 70, 128);
        coexist_channel(15000000, 10, 10, 30, 30, 5, 10, 1500, 310, 70, 128);

    case{9}%varying #device
        zw_time=clock();
        outname='logs/math_No3.txt';

        diary(outname);
        coexist_channel(20000000, 15, 15, 5,  5,  5, 10, 1500, 310, 70, 128);
        coexist_channel(20000000, 15, 15, 10, 10, 5, 10, 1500, 310, 70, 128);
        coexist_channel(20000000, 15, 15, 15, 15, 5, 10, 1500, 310, 70, 128);
        coexist_channel(20000000, 15, 15, 20, 20, 5, 10, 1500, 310, 70, 128);
        coexist_channel(20000000, 15, 15, 25, 25, 5, 10, 1500, 310, 70, 128);
        coexist_channel(20000000, 15, 15, 30, 30, 5, 10, 1500, 310, 70, 128);
    case{10}%varying #device
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        for Packet_w=500:500:1500%[500,1000,1500]
         for aMinBE=[16,32,64]
         for aMaxBE=[256,512,1024]
         for Packet_b=128:-20:48%[128,108,88,68,48]
         for Wc_b=80:-20:40%[70,50,30]
         for Wi_b=320:-80:80%[310,230,150,70]
                    coexist_channel(1000000, wifino, wifino, bmacno, bmacno, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
         end
         end
         end
         end
         end
        end        
        
        diary off;
        
    case{11}%varying everything, for TWC paper.
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        for Packet_w=500:500:1500%[500,1000,1500]
         for aMinBE=[16,32,64]
         for aMaxBE=[256,512,1024]
         for Packet_b=128:-20:48%[128,108,88,68,48]
         for Wc_b=80:-20:40%[70,50,30]
         for Wi_b=320:-80:80%[310,230,150,70]
            for wifino=10:10:30
%                 for bmacno=10:10:80
                    bmacno = 2*wifino;
                    coexist_channel(1000000, wifino, wifino, bmacno, bmacno, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
%                 end
            end
            pause;
         end
         end
         end
         end
         end
        end        
        
        diary off;
    case{12}%varying everything, for TWC paper tuning.
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_qos_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        for Packet_w=1500:500:1500%[500,1000,1500]
         for aMinBE=32:32 % unknown number
         for aMaxBE=1024:1024
         for Packet_b=128:-20:128%[128,108,88,68,48]
         for Wc_b=80:-20:80%[70,50,30]% unknown number
         for Wi_b=160:-80:160%[310,230,150,70]
             for phi = [0.1, 1, 10]
%                 tuning_approximate(1000000, 5, 5, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 10, 10, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 20, 20, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 5, 5, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 10, 10, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 10, 10, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning_approximate(1000000, 20, 20, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 5, 5, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 10, 10, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 20, 20, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 5, 5, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 10, 10, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 10, 10, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
%                 tuning(1000000, 20, 20, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b, phi);
                coexist_channel(1000000, 5, 5, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 10, 10, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 20, 20, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 5, 5, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 10, 10, 20, 20, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 10, 10, 5, 5, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                coexist_channel(1000000, 20, 20, 10, 10, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);

             end
         end
         end
         end
         end
         end
        end        
        
        diary off;
        
	case{13}%varying everything, for Infocom15 paper (saturated case).
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        for Packet_w=1000:500:1000%[500,1000,1500]
         for aMinBE=16%[16,32,64]
         for aMaxBE=[1024] %#ok<*NBRAK>
         for Packet_b=88:-20:88%[128,108,88,68,48]
         for Wc_b=70:-20:70%[70,50,30]
          for Wi_b=310:-80:310%[310,230,150,70]
            for wifino=60:10:60
               for bmacno=40:10:40
                    coexist_channel(1000000, wifino, wifino, bmacno, bmacno, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
                    return;
%                     disp([num2str(Packet_w), ',' , num2str(aMinBE), ',' , num2str(aMaxBE), ',' , ...
%                         num2str(Packet_b), ',' , num2str(Wc_b), ',' , num2str(Wi_b), ',' , ...
%                         num2str(wifino), ',' , num2str(bmacno) ]);
               end
            end
%             pause;
          end
         end
         end
          end
         end
        end        
        
        diary off;
    case{14}%varying everything, for verifying accuracy of the modified way. It seems that it is not accurate, at least when using ns3. 
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        for Packet_w=1000:500:1000%[500,1000,1500]
         for aMinBE=[16,32,64]
         for aMaxBE=[1024] %#ok<*NBRAK>
         for Packet_b=88:-20:88%[128,108,88,68,48]
         for Wc_b=70:-20:30%[70,50,30]
          for Wi_b=310:-80:310%[310,230,150,70]
%             for wifino=20:10:20
%                for bmacno=40:10:40
                    coexist_channel_new(1000000, wifino, wifino, bmacno, bmacno, aMinBE, aMaxBE, Packet_w, Wi_b, Wc_b, Packet_b);
%                     disp([num2str(Packet_w), ',' , num2str(aMinBE), ',' , num2str(aMaxBE), ',' , ...
%                         num2str(Packet_b), ',' , num2str(Wc_b), ',' , num2str(Wi_b), ',' , ...
%                         num2str(wifino), ',' , num2str(bmacno) ]);
%                end
%             end
%             pause;
          end
         end
         end
          end
         end
        end        
        
        diary off;
        
    case{15}%varying everything, for Infocom15 paper (unsaturated case).
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        
        arr_b = 0;arr_w = 0;
        
        for Packet_w=1500:500:1500%[500,1000,1500]
         for Wc_b=70%[30 50 70]%
         for aMaxBE=[1024] %#ok<*NBRAK>
         for Packet_b=48:-20:48%[128,108,88,68,48]
         for aMinBE=16%[16 32 64]%
          for Wi_b=310:-80:310%[310,230,150,70]
            for bmacno=20%[10 20 40 80]%
               for wifino=10%[5 10 20 40]%
%                    for inter_time_b=[500000000:-50000000:50000000]
%                        for inter_time_w=[10000000:-100000:100000]
                            coexist_channel_unsaturated(1000000, wifino, wifino, bmacno, bmacno, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b);
%                        end
%                    end
%                     disp([num2str(Packet_w), ',' , num2str(aMinBE), ',' , num2str(aMaxBE), ',' , ...
%                         num2str(Packet_b), ',' , num2str(Wc_b), ',' , num2str(Wi_b), ',' , ...
%                         num2str(wifino), ',' , num2str(bmacno) ]);
               end
            end
%             pause;
          end
         end
         end
          end
         end
        end        
        
        diary off;

    case{16}%varying everything, for Infocom15 paper tuning for qos.
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

%         diary(outname);
%         diary on;
        
%         arr_b = 0;arr_w = 0;
        
        for Packet_w=1500:500:1500%[500,1000,1500]
         for aMinBE=16%[16 32 64]%
         for aMaxBE=[1024] %#ok<*NBRAK>
         for Packet_b=48:-20:48%[128,108,88,68,48]
         for Wc_b=70%[30 50 70]%
          for Wi_b=310:-80:310%[310,230,150,70]
            for phi = 0.2%[0.2 0.5 1 2 5]%
            for bmacno=20%[5 10 20]%10%
               for wifino=[5 10 20]%5%20%
                for arr_b=1%[1 2 3 4]%1%0.1:0.1:10%10 %per node arrival rate
                    for arr_w=10%[1 10 20 30]%1:1:100 %per node arrival rate%0.01%
                            qos_tuning_unsaturated(1000000, wifino, wifino, bmacno, bmacno, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b, phi);
                    end
                end
%                     disp([num2str(Packet_w), ',' , num2str(aMinBE), ',' , num2str(aMaxBE), ',' , ...
%                         num2str(Packet_b), ',' , num2str(Wc_b), ',' , num2str(Wi_b), ',' , ...
%                         num2str(wifino), ',' , num2str(bmacno) ]);
               end
            end
           end
%             pause;
          end
         end
         end
          end
         end
        end        
        
        diary off;

    case{17}%varying everything, for Infocom15 paper tuning for delay.
        zw_time=clock();
        outname = '';
        outname = sprintf('logs/math_all_%s-%s-%s-%s-%s-%s_%d_%d.txt', num2str(zw_time(2)), num2str(zw_time(3)), ...
            num2str(ceil(zw_time(1))),num2str(zw_time(4)), ...
             num2str(zw_time(5)), num2str(zw_time(6)), wifino, bmacno);

        diary(outname);
        diary on;
        
%         arr_b = 0;arr_w = 0;
        
        for Packet_w=1500:500:1500%[500,1000,1500]
         for aMinBE=16%[16 32 64]%
         for aMaxBE=[1024] %#ok<*NBRAK>
         for Packet_b=48:-20:48%[128,108,88,68,48]
         for Wc_b=70%[30 50 70]%
          for Wi_b=310:-80:310%[310,230,150,70]
            for bmacno=20%[5 10 20]%[10 20 40 80]%20%
               for wifino=[5 10 20]%[5 10 20 40]%10%20%
                for arr_b=1%[1 2 3 4]%1%0.1:0.1:10%10 %per node arrival rate
                    for arr_w=10%[1 10 20 30]%1:1:100 %per node arrival rate%0.01%
                        delay_tuning_unsaturated(1000000, wifino, wifino, bmacno, bmacno, arr_w, aMinBE, aMaxBE, Packet_w, arr_b, Wi_b, Wc_b, Packet_b, 10);
                    end
                end
%                     disp([num2str(Packet_w), ',' , num2str(aMinBE), ',' , num2str(aMaxBE), ',' , ...
%                         num2str(Packet_b), ',' , num2str(Wc_b), ',' , num2str(Wi_b), ',' , ...
%                         num2str(wifino), ',' , num2str(bmacno) ]);
               end
            end
%             pause;
          end
         end
         end
          end
         end
        end        
        
        diary off;        
        
end

diary off;

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

% coexist_channel(500000, 2, 2, 4, 4, 4, 10, 1500, 310, 70, 128);
% coexist_channel(500000, 2, 2, 4, 4, 5, 10, 1500, 310, 70, 128);
% coexist_channel(500000, 2, 2, 4, 4, 6, 10, 1500, 310, 70, 128);
% coexist_channel(500000, 2, 2, 4, 4, 5, 10, 1500, 310, 70, 128);
% coexist_channel(500000, 2, 2, 4, 4, 6, 10, 1500, 310, 70, 128);
ret = 0;

