

function [ out ] = lpl(sim, nodes, timeout, cycle, ratio, cw, traffic )



% clear all;
% clear all;
% clear all;


N = nodes;
%odd id for sender
%even id for receiver.
%receiver_id = sender_id + 1

% t = 0;
L_bmac = 30;%10;%100;%bytes
Ack_timeout_bamc = timeout;%256;
timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot


% mac_queue(1:N) = 0; %# packets
% QUEUE_SIZE = 100;

T = cycle;%ms
rho = ratio;
T_l = T*rho;
% T_s = T-T_l;

timeslot = 10;%us
tic;
sim_time = sim;%ms
true_time = sim_time*1e3/timeslot;

CW(1:2) = [cw 70];
% CW(1:2) = [70 70];

traffic_rate_bmac = traffic;% X = X packets/s


% zw_printf('=============================================\n');
% zw_printf('Simulation time: %g s\n', sim_time/1000);
% zw_printf('# Nodes: %d (%d senders)\n', N, N/2);
% zw_printf('CW: %d time slots\n', CW(1));
% zw_printf('Cycle: %d ms, duty-cycle ratio: %g\n', T, rho);
% zw_printf('Traffic: %d p/s\n', traffic_rate_bmac);
% zw_printf('=============================================\n');
% zw_printf('====ATIM window related:====\n');
% zw_printf('packet size: %d, window size: %d\n', LATIM, ATIM_WINDOW);
% zw_printf('====DATA window related:====\n');
% zw_printf('packet size: %d, DATA window: %d\n', LDATA, DATA_WINDOW);






global nearest_time nearest_time_last rand_save;
global DEBUG_RNG;

nearest_time = 0;
nearest_time_last = -1;
rand_save = -1;

DEBUG_RNG = 0; %debug mode will force the RNG to generate same value for different nodes at a time instance

% global DEBUG_BACKOFF DEBUG_SENSE DEBUG_SENT DEBUG_RECVED DEBUG_WAKEUP DEBUG_SLEEP DEBUG_DR DEBUG_CHANNEL;
global DEBUG_PRT DEBUG_LOG;
DEBUG_PRT = 1;
DEBUG_LOG = 1; % log to file

% DEBUG_BACKOFF = 1;
% DEBUG_SENSE = 1;
% DEBUG_SENT = 1;
% DEBUG_RECVED = 1;
% DEBUG_WAKEUP = 1;
% DEBUG_SLEEP = 1;
% DEBUG_DROP = 1;
% DEBUG_DR = 1;
% DEBUG_CHANNEL = 1;

DEBUG_BACKOFF = 0;
DEBUG_SENSE = 0;
DEBUG_SENT = 0;
DEBUG_RECVED = 0;
DEBUG_WAKEUP = 0;
DEBUG_SLEEP = 0;
DEBUG_DROP = 1;
DEBUG_DR = 0;
DEBUG_CHANNEL = 0;

% DEBUG_SKIP = 1;

SEND_ONE_PACKET = 0;


if DEBUG_LOG

zw_time=clock();
outname=sprintf('sim-%d-%d-%d-%d-%d-%g.txt', ...
    zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
diary(outname);
diary off;
end



rng(1);






channel.idle = 1;
channel.collide = 0;
% zw_printf('channel is idle\n');



% 
% for i=1:100
%     ((-1/traffic_rate_bmac)*log(myrand)*1)
% end
% 
% 
% return;

WAKE_UP = 1;
GOTO_SLEEP = 2;

% DUTY_CYCLE = 1;
% DUTY_CYCLE = 1;
DATA_ARRIVAL = 3;
THE_NEXT = 4;


EVENT_WAITING = 0;%
EVENT_ARRIVAL = 1;%a packet arrives, start backingoff
EVENT_BACKOFF = 2;%starts to backoff
EVENT_SENSE = 3;%backoff done, start carrier sensing 
%(if fail, start another backoff, thus next event is still EVENT_SENSE, o/w, start to send packet, and next event is EVENT_SENT)
EVENT_SENT = 4;%packet has been sent



for ii=1:N
    if mod(ii,2) % transmitter
        node_struct(ii).waiting = 1; 
        node_struct(ii).next_event = -1;
        node_struct(ii).next_event_time = Inf;
        node_struct(ii).wakeup_event_time = Inf;
        node_struct(ii).sleep_event_time = Inf;
        
        node_struct(ii).dtar_event_time = ceil((-1/traffic_rate_bmac)*log(myrand)*1e5);%s to 10 us;
        zw_printf('[%d @ %d us] wait packet to arrive at %d\n', ii, nearest_time*10, node_struct(ii).dtar_event_time);
        node_struct(ii).cw_stage = 1;
%         node_struct(ii).collide = 0;
        node_struct(ii).awake = 0;
        node_struct(ii).awaketime = 0;
        node_struct(ii).last_awake = 0;
        node_struct(ii).recv_awake = 0;
        node_struct(ii).queue = 0;
        node_struct(ii).success = 0;
        node_struct(ii).collision = 0; %fail due to collision
        node_struct(ii).sleep = 0; %fail due to sleeping
        node_struct(ii).drop = 0;
        node_struct(ii).count = 0;
        node_struct(ii).txtime = 0;
        node_struct(ii).txstart = 0;
    else % receiver of ii-1
        node_struct(ii).waiting = 0; 
        node_struct(ii).next_event = -1;
        node_struct(ii).next_event_time = Inf;
%         node_struct(ii).ducy_event_time = int32(T/2*100*myrand);
        node_struct(ii).wakeup_event_time = int32(T/2*100*myrand);
        node_struct(ii).sleep_event_time = Inf;
        node_struct(ii).dtar_event_time = Inf;
        zw_printf('[%d @ %d us] wait node to wakeup at %d\n', ii, nearest_time*10, node_struct(ii).wakeup_event_time);
        node_struct(ii).cw_stage = -1;
%         node_struct(ii).collide = 0;
        node_struct(ii).awake = 0;
        node_struct(ii).awaketime = 0;
        node_struct(ii).last_awake = 0;
        node_struct(ii).recv_awake = -1;
        node_struct(ii).queue = -1;
        node_struct(ii).success = -1;
        node_struct(ii).collision = -1;
        node_struct(ii).sleep = -1; 
        node_struct(ii).drop = -1;
        node_struct(ii).count = 0;
    end
end

%     node_struct(4).next_event = EVENT_WAITING;
%     node_struct(4).next_event_time = 1;
%     node_struct(4).dtar_event_time = 0;
%     node_struct(3).next_event = EVENT_ARRIVAL;
%     node_struct(3).next_event_time = 1;
%     node_struct(3).ducy_event_time = 0;
%     node_struct(2).next_event = EVENT_SENT;
%     node_struct(2).next_event_time = 0;
%     node_struct(1).next_event = EVENT_SENSE;
%     node_struct(1).next_event_time = 0;
%     node_struct(1).ducy_event_time = 0;



% next_event(1:N) = EVENT_WAITING;
% % next_event(4) = EVENT_WAITING;
% % next_event(3) = EVENT_ARRIVAL;
% % next_event(2) = EVENT_SENSE;
% % next_event(1) = EVENT_SENT;
% 
% 
% next_event_time(1:N) = 0;
% % next_event_time(3:N) = 1;
% 
% ducy_event_time(1:N) = Inf;%duty cycling event time
% dtar_event_time(1:N) = Inf;%data arrival event time









while 1
% for ii=1:50
    
%     next_time = min(next_event_time);
% 
%     ducy_time = min(ducy_event_time);
% 
%     dtar_time = min(dtar_event_time);

    nearest_time = min([node_struct.wakeup_event_time node_struct.sleep_event_time node_struct.dtar_event_time node_struct.next_event_time]); %find the nearest time
    nearest_list = find([node_struct.wakeup_event_time node_struct.sleep_event_time node_struct.dtar_event_time node_struct.next_event_time]==nearest_time); %which one is the nearest
    % 1~N: next event, N+1~2N: duty cycle event, 2N+1, 3N

%     zw_printf('nearest_time: %d\n', nearest_time);
%     zw_printf('nearest list: \n');
%     disp(nearest_list);
    
    
    nearest_type = (ceil(nearest_list/(N))); %check the even type and remove duplicates
    %  1:wakeup event; 2: sleep event; 3: data arrival event; 4: next event;
    
%     zw_printf('nearest_type list: \n');
%     disp(nearest_type);
    
    
    nearest_from = mod(nearest_list-1,N)+1;
%     %this makes sure that EVENT_SENT (=3) lies before
%     %EVENT_SENSE (=2), thus if a node finishes its transmission
%     %at the same time with some nodes who start to sense the
%     %channel, the latter will have consistant behaviour because the channel will be set to idle before any of the nodes sensing the channel.
%     A = nearest_from(nearest_type==THE_NEXT);
%     B = [A; node_struct(A).next_event];
%     [Y,I]=sort(B(2,:), 'descend');
%     nearest_from(nearest_type==THE_NEXT) = B(1,I);

%     zw_printf('nearest_from list: \n');
%     disp(nearest_from);
    
    nearest_info = [nearest_from; nearest_type];
    
%     return;

    
    jj = 0;
    for ii=nearest_info(1,:)
        jj = jj + 1;
    
%         zw_printf('nearest_from: %d\n', ii);

        switch nearest_info(2,jj)

            case THE_NEXT %next event
%                 zw_printf('THE_NEXT event\n');
%                 re = find(next_event_time==nearest_time); %list of nodes who have the nearest event
%                 zw_printf('next event list:\n');
%                 disp(re);
%                 next_event(re)
%                 re = sort(re, 'descend'); 
% %                 return;
%                 
%                 for jj=re
                    switch node_struct(ii).next_event
                        
%                         case EVENT_WAITING
%                             node_struct(ii).next_event_time = node_struct(ii).next_event_time + ceil((-1/traffic_rate_bmac)*log(myrand)*1e5)*timeslot_ratio;%s to 10 us
%                             zw_printf('[%d @ %d] wait packet to arrive at %d\n', ii, nearest_time, node_struct(ii).next_event_time);
%                             node_struct(ii).next_event = EVENT_ARRIVAL;
%                             
%                         case EVENT_ARRIVAL
%                             node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/traffic_rate_bmac)*log(myrand)*1e5)*timeslot_ratio;
%                             node_struct(ii).ducy_event_time = node_struct(ii).ducy_event_time + T*1e2;%ms to 10us
%                             node_struct(ii).next_event_time = node_struct(ii).next_event_time + 1;
%                             node_struct(ii).next_event = EVENT_BACKOFF;
%                             zw_printf('[%d @ %d] packet arrives, backoff will begin at %d\n', ii, nearest_time, node_struct(ii).next_event_time);
% %                             node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(10+myrand*CW(node_struct(ii).cw_stage))*timeslot_ratio;%backoff
% %                             zw_printf('[%d] backingoff done (sensing starts) at %d\n', ii, node_struct(ii).next_event_time);
% %                             node_struct(ii).next_event = EVENT_SENSE;
                            
                        case EVENT_BACKOFF
%                             if node_struct(ii).awake==0
%                                 if DEBUG_SKIP; zw_printf('[%d @ %d us] has been reset, skip backing off\n', ii, nearest_time*10); end
%                                 continue;
%                             end
                            node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(10+myrand*CW(node_struct(ii).cw_stage))*timeslot_ratio;%backoff
                            node_struct(ii).next_event = EVENT_SENSE;
                            if DEBUG_BACKOFF; zw_printf('[%d @ %d us] backoff starts, CCA will start at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end

                        case EVENT_SENSE
                    %         zw_printf('%d\n', i);
                    % check the channel
%                             if node_struct(ii).awake==0
%                                 if DEBUG_SKIP; zw_printf('[%d @ %d us] has been reset, skip CCA\n', ii, nearest_time*10); end
%                                 continue;
%                             end
                            
                            if channel.idle==0%channel is busy

                                % node_struct(ii).cw_stage = min(node_struct(ii).cw_stage+1, 2); %assume there is only one stage
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(10+myrand*CW(node_struct(ii).cw_stage))*timeslot_ratio;%backoff
                                node_struct(ii).next_event = EVENT_SENSE;
                                if DEBUG_SENSE; zw_printf('[%d @ %d us] backoff done, CCA fails, and will restart at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                            else
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(L_bmac)*timeslot_ratio;%backoff
                                node_struct(ii).next_event = EVENT_SENT;
                                if DEBUG_SENSE; zw_printf('[%d @ %d us] backoff done, CCA succeeds, transmitting and will be done at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                % check if the corresponding receiver is awake
                                % if yes, set node_struct.recv_awake = 1.
                                % o/w, clear it to 0.
                                if node_struct(ii+1).awake
                                    node_struct(ii).recv_awake = 1;
                                else
                                    node_struct(ii).recv_awake = 0;
                                end
                            end
                            
                        case EVENT_SENT
                            
                            % see if collided
                            % if yes
                            %       clear node_struct.recv_awake
                            % o/w
                            %       if the receiver was awake (recorded)
                            %           check again to see if it is still awake
                            %           if yes, the transmission succeeds
                            %           o/w, the transmission fails.
                            %       o/w, the transmission fails.
                            
                            if channel.collide==1 % see if collided
                                % the transmission fails.
                                % clear node_struct.recv_awake
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ack_timeout_bamc)*timeslot_ratio;%backoff
                                node_struct(ii).next_event = EVENT_BACKOFF;
                                if DEBUG_SENT; zw_printf('[%d @ %d us] sending fails (collision), wait until Ack timesout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                node_struct(ii).collision = node_struct(ii).collision + 1;
                            else
                                if node_struct(ii).recv_awake         % if the receiver was awake (recorded)

                                    if node_struct(ii+1).awake        % check again to see if it is still awake, if yes, the transmission succeeds
                                        
%                                         zw_printf('[%d @ %d] sending succeed, sleep immediately \n', ii, nearest_time);
%                                         node_struct(ii).sleep_event_time = Inf;%ms to 10us
%                                         node_struct(ii).next_event = -1;
%                                         node_struct(ii).next_event_time = Inf;
%                                         node_struct(ii).awake = 0;
                                        
                                        node_struct(ii).queue = node_struct(ii).queue - 1;
                                        node_struct(ii).count = node_struct(ii).count + 1;
                                        node_struct(ii).success = node_struct(ii).success + 1;
                                        
                                        if node_struct(ii).queue ~= 0 %has a packet in the queue
                                            node_struct(ii).sleep_event_time = nearest_time +1+ T*100;%ms to 10us
                                            node_struct(ii).next_event_time = nearest_time + 1;
                                            node_struct(ii).next_event = EVENT_BACKOFF;
                                            node_struct(ii).waiting = 0;
                                            node_struct(ii).awake = 1;
                                            node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                                            node_struct(ii).txstart = nearest_time;
                                            if DEBUG_SENT; zw_printf('[%d @ %d us] sending succeed, and there is a new packet, begin to backoff at %d us and will fall sleep at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10, node_struct(ii).sleep_event_time*10); end
                                        else
                                            node_struct(ii).sleep_event_time = Inf;%ms to 10us
                                            node_struct(ii).next_event = -1;
                                            node_struct(ii).next_event_time = Inf;
                                            node_struct(ii).awake = 0;
                                            node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                            node_struct(ii).waiting = 1;
                                            node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                                            
                                            if DEBUG_SENT; zw_printf('[%d @ %d us] sending succeed, and there is no new packet, sleep immediately\n', ii, nearest_time*10); end
                                        end
                        
                                        node_struct(ii+1).awake = 0;
                                        node_struct(ii+1).sleep_event_time = Inf;
                                        node_struct(ii+1).wakeup_event_time = nearest_time + (T-T_l)*100; %important fix!!!! 6/22/2016
                                        if DEBUG_RECVED; zw_printf('[%d @ %d us] recved succeed, sleep immediately and will wakeup at %d us\n', ii+1, nearest_time*10, node_struct(ii+1).wakeup_event_time*10); end
                                        
                                    else        % o/w, the transmission fails.
                                        node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ack_timeout_bamc)*timeslot_ratio;%backoff
                                        node_struct(ii).next_event = EVENT_BACKOFF;
                                        if DEBUG_SENT; zw_printf('[%d @ %d us] sending fails (awake to sleep), wait until Ack timesout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                        node_struct(ii).sleep = node_struct(ii).sleep + 1;
                                    end
                                else          % o/w, the transmission fails.
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ack_timeout_bamc)*timeslot_ratio;%backoff
                                    node_struct(ii).next_event = EVENT_BACKOFF;
                                    if DEBUG_SENT; zw_printf('[%d @ %d us] sending fails (sleeping), wait until Ack timesout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                    node_struct(ii).sleep = node_struct(ii).sleep + 1;
                                    
                                    
                                end
                                
                            end
                            
                            
                    end

%                 end

            case WAKE_UP % wakeup event
%                 zw_printf('WAKEUP event\n');
                
                if mod(ii,2) % transmitter
                    error(sprintf('[%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10));
%                     if node_struct(ii).next_event == EVENT_SENT % is currently sending
%                         node_struct(ii).sleep_event_time = node_struct(ii).next_event_time + 1;%ms to 10us
%                         zw_printf('[%d @ %d] node is sending, defer to sleep at %d\n', ii, nearest_time, node_struct(ii).sleep_event_time);
%                         
%                     else
%                         node_struct(ii).sleep_event_time = Inf;%ms to 10us
%                         node_struct(ii).next_event = -1;
%                         node_struct(ii).next_event_time = Inf;
%                         zw_printf('[%d @ %d] node is not sending, sleep immediately\n', ii, nearest_time);
%                         node_struct(ii).awake = 0;
%                     end
                    
                else % receiver
                    if node_struct(ii).awake == 0
                        node_struct(ii).wakeup_event_time = node_struct(ii).wakeup_event_time + T*100;%ms to 10us
                        if node_struct(ii).sleep_event_time == Inf; node_struct(ii).sleep_event_time = nearest_time; end
                        node_struct(ii).sleep_event_time = node_struct(ii).sleep_event_time + T_l*100;%ms to 10us
                        if DEBUG_WAKEUP; zw_printf('[%d @ %d us] wakes up, will go to sleep at %d us\n', ii, nearest_time*10, node_struct(ii).sleep_event_time*10); end
                        node_struct(ii).awake = 1;
                    else
                        
                        
                        if DEBUG_WAKEUP; error(sprintf('[%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10)); end
%                         node_struct(ii).ducy_event_time = node_struct(ii).ducy_event_time + T*10;%*1e2;%ms to 10us
%                         
%                         zw_printf('[%d @ %d] goes to sleep, will wake up at %d\n', ii, nearest_time, node_struct(ii).ducy_event_time);
%                         node_struct(ii).awake = 0;
                       
                    end
                end
                
            case GOTO_SLEEP % sleep event
%                 zw_printf('GOTO SLEEP event\n');
                
                if mod(ii,2) % transmitter
                    
                    if node_struct(ii).next_event == EVENT_SENT % is currently sending
                        node_struct(ii).sleep_event_time = node_struct(ii).next_event_time + 1;%ms to 10us
                        if DEBUG_SLEEP; zw_printf('[%d @ %d us] node is sending, defer to sleep at %d us\n', ii, nearest_time*10, node_struct(ii).sleep_event_time*10); end
                        
                    else
                        node_struct(ii).queue = node_struct(ii).queue - 1;
                        node_struct(ii).count = node_struct(ii).count + 1;
                        node_struct(ii).drop = node_struct(ii).drop + 1;
                        
                        if node_struct(ii).queue ~= 0 %has a packet in the queue
                            node_struct(ii).sleep_event_time = nearest_time +1+ T*100;%ms to 10us
                            node_struct(ii).next_event_time = nearest_time + 1;
                            node_struct(ii).next_event = EVENT_BACKOFF;
                            node_struct(ii).waiting = 0;
                            node_struct(ii).awake = 1;
                            node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                            node_struct(ii).txstart = nearest_time;
                            if DEBUG_DROP; zw_printf('[%d @ %d us] node is not sending, and there is a new packet, drop and begin to backoff at %d us and will fall sleep at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10, node_struct(ii).sleep_event_time*10); end
                        else
                            node_struct(ii).sleep_event_time = Inf;%ms to 10us
                            node_struct(ii).next_event = -1;
                            node_struct(ii).next_event_time = Inf;
                            if DEBUG_DROP; zw_printf('[%d @ %d us] node is not sending, and there is no new packet, drop and sleep immediately\n', ii, nearest_time*10); end
                            node_struct(ii).awake = 0;
                            node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                            node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                            node_struct(ii).waiting = 1;
                        end
                    end
                    
                else % receiver
                    if node_struct(ii).awake == 0
                        error(sprintf('[%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10));
%                         node_struct(ii).ducy_event_time = node_struct(ii).ducy_event_time + T*10;%*1e2;%ms to 10us
%                         zw_printf('[%d @ %d] wakes up, will go to sleep at %d\n', ii, nearest_time, node_struct(ii).ducy_event_time);
%                         node_struct(ii).awake = 1;
                    else
                        node_struct(ii).sleep_event_time = Inf;
                        if DEBUG_SLEEP; zw_printf('[%d @ %d us] goes to sleep without receiving successfully\n', ii, nearest_time*10); end
                        
%                         if DEBUG_SLEEP; zw_printf('[%d @ %d us] goes to sleep, will wake up at %d us\n', ii, nearest_time*10, node_struct(ii).wakeup_event_time*10); end
                        node_struct(ii).awake = 0;
                       
                    end
                end
                
                
                

            case DATA_ARRIVAL %data arrival event
%                 zw_printf('DATA_ARRIVAL event\n');
                if DEBUG_DR; zw_printf('[%d @ %d us] packet arrives..., ', ii, nearest_time*10); end
                if node_struct(ii).waiting == 1
                    node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/traffic_rate_bmac)*log(myrand)*1e5);
                    if SEND_ONE_PACKET; node_struct(ii).dtar_event_time = Inf; end
                    if DEBUG_DR; zw_printf('next packet will arrive at %d us, ', node_struct(ii).dtar_event_time*10); end
                    if node_struct(ii).sleep_event_time == Inf; node_struct(ii).sleep_event_time = nearest_time; end
                    node_struct(ii).sleep_event_time = node_struct(ii).sleep_event_time +1 + T*100;%ms to 10us
                    if DEBUG_DR; zw_printf('node will sleep at %d us, ', node_struct(ii).sleep_event_time*10); end
                    if node_struct(ii).next_event_time == Inf; node_struct(ii).next_event_time = nearest_time; end
                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + 1;
                    if DEBUG_DR; zw_printf('backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
                    node_struct(ii).next_event = EVENT_BACKOFF;
                    node_struct(ii).waiting = 0;
                    node_struct(ii).awake = 1;
                    node_struct(ii).last_awake = nearest_time;
                    node_struct(ii).queue = node_struct(ii).queue + 1;
                    node_struct(ii).txstart = nearest_time;
                else
                    node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/traffic_rate_bmac)*log(myrand)*1e5);
                    if SEND_ONE_PACKET; node_struct(ii).dtar_event_time = Inf; end
                    if DEBUG_DR; zw_printf('node is sending, put the packet into the queue\n'); end
                    node_struct(ii).queue = node_struct(ii).queue + 1;
                end

        end    
    
    
    end
    
    transmittings = find([node_struct.next_event] == EVENT_SENT);
%     disp(transmittings);
    if isempty(transmittings)
        channel.idle = 1;
        if DEBUG_CHANNEL; zw_printf('channel is idle\n'); end
    else
        channel.idle = 0;
        if DEBUG_CHANNEL; zw_printf('channel is busy...'); end
        if size(transmittings,2) > 1 %multiple transmitings.
            if DEBUG_CHANNEL; zw_printf('and there is a collision\n'); end
            channel.collide = 1;
        else
            if DEBUG_CHANNEL; zw_printf('and transmission is collision free\n'); end
            channel.collide = 0;
        end
        
    end
    
    
    
    
    t = nearest_time;
    
    if t>=true_time
        
       break;
       
    end
    
    
end

toc;
if DEBUG_LOG
    diary off;
end

outname=sprintf('stat(%d-%d-%d-%g-%d-%g).mat', sim_time/1000, N, T, rho, CW(1), traffic_rate_bmac);

save(outname, 'node_struct');

out = node_struct;

end
