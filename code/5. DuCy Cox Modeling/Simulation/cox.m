

function out = cox(sim, nodes_B, packet_B, timeout_B, cycle_B, ratio_B, cw_B, traffic_B, ...
                            nodes_W, packet_W, timeout_W, cycle_W, ratio_W, cw_W, traffic_W)
global showdetailed_nt;
global afig singlefig;
global asfig;

global col;

global packet_mactime;
global packet_sucinatim;


% clear all;
% clear all;
% clear all;

timeslot = 10;

N_W = nodes_W;
LA_W = 1;
L_W = packet_W;%*8/54/timeslot;%1000 bytes to WiFi time slots = 15.
Ackto_W = timeout_W;%256; Ack timeout
T_W = cycle_W;%ms
Rho_W = ratio_W;
CW_W(1:2) = [cw_W 16];
TR_W = traffic_W;% X = X packets/s

Actnode_cur(1:T_W*1e3/timeslot)=-1;
Actnode_cur(1)=N_W/2;
last_offset = -1; % relative position within T_W
last_position = -1; % absolute position within T_W

Actnode_stat(1:T_W*1e3/timeslot)=0;

N_B = nodes_B;
%odd id for sender
%even id for receiver.
%receiver_id = sender_id + 1

% t = 0;
Ackto_B = timeout_B;%256; Ack timeout
timeslot_ratio = 3;%bmac timeslot is 3 x wifi timeslot
L_B = packet_B;%*8/0.25/timeslot_ratio/timeslot;%94bytes to BMAC time slots = 100 (i.e 300 WiFi time slots).


% mac_queue(1:N) = 0; %# packets
% QUEUE_SIZE = 100;

T_B = cycle_B;%ms
Rho_B = ratio_B;
Tl_B = T_B*Rho_B; %listen time
% T_s = T-T_l;

tic;
sim_time = sim;%ms
true_time = sim_time*1e3/timeslot;

CW_B(1:2) = [cw_B 70];
% CW(1:2) = [70 70];

TR_B = traffic_B;% X = X packets/s


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
% DEBUG_PRT = 1;
DEBUG_LOG = 0; % log to file

DEBUG_DETAIL = 0;

if DEBUG_DETAIL

    DEBUG_BACKOFF = 1;
    DEBUG_SENSE = 1;
    DEBUG_SENT = 1;
    DEBUG_COLL = 1;
    DEBUG_RECVED = 1;
    DEBUG_WAKEUP = 1;
    DEBUG_SLEEP = 1;
    DEBUG_DR = 1;
    DEBUG_CHANNEL = 1;
else

    DEBUG_BACKOFF = 0;
    DEBUG_SENSE = 0;
    DEBUG_SENT = 0;
    DEBUG_COLL = 0;
    DEBUG_RECVED = 0;
    DEBUG_WAKEUP = 0;
    DEBUG_SLEEP = 0;
    DEBUG_DR = 0;
    DEBUG_CHANNEL = 0;

end

% DEBUG_SKIP = 1;

SEND_ONE_PACKET = 0;
SATURATE_LPL = 0;
SATURATE_PSM = 0;

if DEBUG_LOG

    zw_time=clock();
    outname=sprintf('sim-%d-%d-%d-%d-%d-%g.txt', ...
        zw_time(2), zw_time(3), ceil(zw_time(1)), zw_time(4), zw_time(5), zw_time(6));
    diary(outname);
    diary off;
end

rng(1);




% 
% for i=1:100
%     ((-1/traffic_rate_bmac)*log(myrand)*1)
% end
% 
% 
% return;

WAKE_UP = 1;
GOTO_SLEEP = 2;
DATA_ARRIVAL = 3;
THE_NEXT = 4;


% EVENT_WAITING = 0;%
% EVENT_ARRIVAL = 1;%a packet arrives, start backingoff
EVENT_NONE = -1;
EVENT_BACKOFF = 2;%starts to backoff (DIFS is done)
EVENT_SENSE = 3;%backoff done, start carrier sensing 
%(if fail, start another backoff, thus next event is still EVENT_SENSE, o/w, start to send packet, and next event is EVENT_SENT)
EVENT_SENT = 4;%packet has been sent


DEVICE_LPL = 1;
DEVICE_PSM = 2;

DIFS = 3;

% initial_ATIM_starting = myrand * T_W;
msto10us = 1e2; %ms to 10us

% disp(initial_ATIM_starting);



channel.idle = 1;
channel.collide = 0;
% channel.idle2busy = 0; % when channel becomes busy, this becomes 1 at that time instance, and then it becomes 0.
% channel.busyfor = 0;
channel.busyuntil = Inf;
% zw_printf('channel is idle\n');
channel.atimstart = -1;
channel.datastart = -1;



for ii=1:N_W/2 %only considering the transmitters for PSM
    node_struct(ii).type = DEVICE_PSM; 
    node_struct(ii).waiting = 1; 
    node_struct(ii).next_event = EVENT_NONE;
    node_struct(ii).next_event_time = Inf;
    node_struct(ii).wakeup_event_time = 0;%ceil(initial_ATIM_starting*msto10us); %ATIM window starts (DATA window ends)
    node_struct(ii).sleep_event_time = Inf;%initial_ATIM_starting + Rho_W*T_W*msto10us; %ATIM window ends (DATA window starts)

    if SATURATE_PSM == 1
        node_struct(ii).sat = 1;
        node_struct(ii).dtar_event_time = Inf;
        zw_printf('[PSM:%d @ %d us] saturated traffic!\n', ii, nearest_time*10);
    else
        node_struct(ii).sat = 0;
        node_struct(ii).dtar_event_time = ceil((-1/TR_W)*log(myrand)*1e5);%s to 10 us;
        zw_printf('[PSM:%d @ %d us] wait packet to arrive at %d\n', ii, nearest_time*10, node_struct(ii).dtar_event_time*10);
    end
    node_struct(ii).cw_stage = 0;
%         node_struct(ii).collide = 0;
    node_struct(ii).awake = 0;
    node_struct(ii).awaketime = 0;
    node_struct(ii).last_awake = 0;
    node_struct(ii).recv_awake = 0;
    node_struct(ii).queue = 0;
    node_struct(ii).success = 0; %statistics for number of successes.
    node_struct(ii).collision = 0; %statistics for number of collisions.
    node_struct(ii).sleep = 0; %statistics for number of sleeps.
    node_struct(ii).drop = 0;
    node_struct(ii).count = 0;
    node_struct(ii).txtime = 0;
    node_struct(ii).txstart = 0;
%     node_struct(ii).difsing = 0;
%     node_struct(ii).backingoff = 0; %backoff flag
    node_struct(ii).sleepindw = 0; %sleep in data window when no data in atim
    node_struct(ii).backoffleft = -1;%number of backoff left (when freeze happens)
    node_struct(ii).succeeded = 0; % succeed in atim or data window
    node_struct(ii).inatim = 1; % whether in atim 
    node_struct(ii).waitingchannel = 1;
    
    
    
end

TRANSMITTER = 1;
RECEIVER = 2;

kk=0;
for ii=N_W/2+1:N_W/2+N_B
    node_struct(ii).type = DEVICE_LPL; 
    kk = kk + 1;
    if mod(kk,2) % transmitter
        node_struct(ii).device = TRANSMITTER; 
        node_struct(ii).waiting = 1; 
        node_struct(ii).next_event = EVENT_NONE;
        node_struct(ii).next_event_time = Inf;
        node_struct(ii).wakeup_event_time = Inf;
        node_struct(ii).sleep_event_time = Inf;
        
        if SATURATE_LPL == 1
            node_struct(ii).sat = 1;
            node_struct(ii).dtar_event_time = 0;%myrand*T_B*1e1;
            zw_printf('[LPL:%d @ %d us] saturated traffic!\n', ii, nearest_time*10);
        else
            node_struct(ii).sat = 0;
            node_struct(ii).dtar_event_time = ceil((-1/TR_B)*log(myrand)*1e5);%s to timeslot;
            zw_printf('[LPL:%d @ %d us] wait packet to arrive at %d\n', ii, nearest_time*10, node_struct(ii).dtar_event_time*10);
        end
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
        node_struct(ii).txtime = 0; %total time 
        node_struct(ii).txstart = 0; %time 
        
        
        
    else % receiver of ii-1
        node_struct(ii).device = RECEIVER; 
        node_struct(ii).waiting = 0; 
        node_struct(ii).next_event = EVENT_NONE;
        node_struct(ii).next_event_time = Inf;
%         node_struct(ii).ducy_event_time = int32(T/2*100*myrand);
        node_struct(ii).wakeup_event_time = int32(T_B/2*100*myrand);
        node_struct(ii).sleep_event_time = Inf;
        node_struct(ii).dtar_event_time = Inf;
        zw_printf('[LPL:%d @ %d us] wait node to wakeup at %d\n', ii, nearest_time*10, node_struct(ii).wakeup_event_time*10);
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


N=N_W/2+N_B;




Round = 0;
current_round = 0;
print_density = 200;

while 1
% for ii=1:50
    
%     next_time = min(next_event_time);
% 
%     ducy_time = min(ducy_event_time);
% 
%     dtar_time = min(dtar_event_time);

    nearest_time = min([node_struct.wakeup_event_time node_struct.sleep_event_time node_struct.dtar_event_time node_struct.next_event_time]); %find the nearest time
    nearest_list = find([node_struct.wakeup_event_time node_struct.sleep_event_time node_struct.dtar_event_time node_struct.next_event_time]==nearest_time); %which one is the nearest
    % 1~N: wakeup event, N+1~2N: sleep event, 2N+1~3N: data arrival event, 3N+1~4N: next event

%     zw_printf('nearest_time: %d\n', nearest_time*10);
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
                            if node_struct(ii).type == DEVICE_LPL
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(10+myrand*CW_B(node_struct(ii).cw_stage))*timeslot_ratio;%backoff
                                node_struct(ii).next_event = EVENT_SENSE;
                                if DEBUG_BACKOFF; zw_printf('[LPL:%d @ %d us] backoff starts, CCA will start at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                            else %PSM
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + node_struct(ii).backoffleft;%int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
                                if DEBUG_BACKOFF; zw_printf('[PSM:%d @ %d us] DIFS is done, backoff starts and will be done at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                node_struct(ii).next_event = EVENT_SENSE; %for PSM, we just reuse EVENT_SENSE because it will send the 
                                                                          %data immediately
%                                 node_struct(ii).difsing = 0;
%                                 node_struct(ii).backingoff = 1;
                                node_struct(ii).backoffleft = -1;
                            end
                        case EVENT_SENSE
                    %         zw_printf('%d\n', i);
                    % check the channel
%                             if node_struct(ii).awake==0
%                                 if DEBUG_SKIP; zw_printf('[%d @ %d us] has been reset, skip CCA\n', ii, nearest_time*10); end
%                                 continue;
%                             end
                            if node_struct(ii).type == DEVICE_LPL
                            
                                if channel.idle==0%channel is busy

                                    % node_struct(ii).cw_stage = min(node_struct(ii).cw_stage+1, 2); %assume there is only one stage
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(10+myrand*CW_B(node_struct(ii).cw_stage))*timeslot_ratio;%backoff
                                    node_struct(ii).next_event = EVENT_SENSE;
                                    if DEBUG_SENSE; zw_printf('[LPL:%d @ %d us] backoff done, CCA fails, and will backoff again and retry at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                else
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(L_B)*timeslot_ratio;%transmission
                                    node_struct(ii).next_event = EVENT_SENT;
%                                     channel.idle2busy = 1;
%                                     channel.busyfor = max(int32(L_B)*timeslot_ratio, channel.busyfor);
%                                     channel.busyuntil = nearest_time + channel.busyfor;
                                    if DEBUG_SENSE; zw_printf('[LPL:%d @ %d us] backoff done, CCA succeeds, transmitting and will be done at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                    % check if the corresponding receiver is awake
                                    % if yes, set node_struct.recv_awake = 1.
                                    % o/w, clear it to 0.
                                    if node_struct(ii+1).awake
                                        node_struct(ii).recv_awake = 1;
                                    else
                                        node_struct(ii).recv_awake = 0;
                                    end
                                end
                            else %PSM, at this point, the PSM will send the packet immediately
                                
                                if node_struct(ii).inatim
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(LA_W);%transmission
%                                     channel.busyfor = max(int32(LA_W), channel.busyfor);
%                                     channel.busyuntil = nearest_time + channel.busyfor;
                                else
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(L_W);%transmission
%                                     channel.busyfor = max(int32(L_W), channel.busyfor);
%                                     channel.busyuntil = nearest_time + channel.busyfor;
                                end
                                node_struct(ii).next_event = EVENT_SENT;
%                                 channel.idle2busy = 1;
                                if DEBUG_SENSE; zw_printf('[PSM:%d @ %d us] backoff done, transmitting and will be done at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
%                                 node_struct(ii).backingoff = 0;
                                
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
                           if node_struct(ii).type == DEVICE_LPL
                            
                                if channel.collide==1 % see if collided
                                    % the transmission fails.
                                    % clear node_struct.recv_awake
                                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ackto_B)*timeslot_ratio;%backoff
                                    node_struct(ii).next_event = EVENT_BACKOFF;
                                    if DEBUG_COLL; zw_printf('[LPL:%d @ %d us] sending fails (collision), wait until Ack timeout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                    node_struct(ii).collision = node_struct(ii).collision + 1;
                                else
                                    if node_struct(ii).recv_awake         % if the receiver was awake (recorded)

                                        if node_struct(ii+1).awake        % check again to see if it is still awake, if yes, the transmission succeeds

    %                                         zw_printf('[%d @ %d] sending succeed, sleep immediately \n', ii, nearest_time);
    %                                         node_struct(ii).sleep_event_time = Inf;%ms to 10us
    %                                         node_struct(ii).next_event = EVENT_NONE;
    %                                         node_struct(ii).next_event_time = Inf;
    %                                         node_struct(ii).awake = 0;

                                            node_struct(ii).queue = node_struct(ii).queue - 1;
                                            node_struct(ii).count = node_struct(ii).count + 1;
                                            node_struct(ii).success = node_struct(ii).success + 1;

%                                             assert(node_struct(ii).queue >= 0);
                                            if node_struct(ii).queue > 0 || node_struct(ii).sat == 1%has a packet in the queue
                                                node_struct(ii).sleep_event_time = nearest_time +1+ T_B*msto10us;%ms to 10us
                                                node_struct(ii).next_event_time = nearest_time + 1;
                                                node_struct(ii).next_event = EVENT_BACKOFF;
                                                node_struct(ii).waiting = 0;
                                                node_struct(ii).awake = 1;
                                                node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                                node_struct(ii).last_awake = nearest_time;
                                                node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                                                node_struct(ii).txstart = nearest_time;
                                                if DEBUG_SENT; zw_printf('[LPL:%d @ %d us] sending succeed, and there is a new packet, begin to backoff at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                            else
                                                node_struct(ii).sleep_event_time = Inf;%ms to 10us
                                                node_struct(ii).next_event = EVENT_NONE;
                                                node_struct(ii).next_event_time = Inf;
                                                node_struct(ii).awake = 0;
                                                node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                                node_struct(ii).waiting = 1;
                                                node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;

                                                if DEBUG_SENT; zw_printf('[LPL:%d @ %d us] sending succeed, and there is no new packet, sleep immediately\n', ii, nearest_time*10); end
                                            end

                                            node_struct(ii+1).awake = 0;
                                            node_struct(ii+1).sleep_event_time = Inf;
                                            node_struct(ii+1).wakeup_event_time = nearest_time + (T_B-Tl_B)*100; %important fix!!!! 6/22/2016
                                            if DEBUG_RECVED; zw_printf('[LPL:%d @ %d us] recved succeed, sleep immediately and will wakeup at %d us\n', ii+1, nearest_time*10, node_struct(ii+1).wakeup_event_time*10); end

                                        else        % o/w, the transmission fails.
                                            node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ackto_B)*timeslot_ratio;%backoff
                                            node_struct(ii).next_event = EVENT_BACKOFF;
                                            if DEBUG_SENT; zw_printf('[LPL:%d @ %d us] sending fails (awake to sleep), wait until Ack timeout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                            node_struct(ii).sleep = node_struct(ii).sleep + 1;
                                        end
                                    else          % o/w, the transmission fails.
                                        node_struct(ii).next_event_time = node_struct(ii).next_event_time + int32(Ackto_B)*timeslot_ratio;%backoff
                                        node_struct(ii).next_event = EVENT_BACKOFF;
                                        if DEBUG_SENT; zw_printf('[LPL:%d @ %d us] sending fails (sleeping), wait until Ack timeout at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                                        node_struct(ii).sleep = node_struct(ii).sleep + 1;


                                    end

                                end
                            
                           else %PSM
                               
                                if channel.collide==1 % see if collided
                                    node_struct(ii).next_event_time = Inf;%backoff
                                    node_struct(ii).next_event = EVENT_NONE;
                                    if DEBUG_COLL; zw_printf('[PSM:%d @ %d us] sending fails (collision), wait until window reset or channel idle, ', ii, nearest_time*10); end
                                    node_struct(ii).collision = node_struct(ii).collision + 1;
%                                     node_struct(ii).difsing = 1;
                                    node_struct(ii).cw_stage = min(node_struct(ii).cw_stage + 1, 5);
                                    if DEBUG_COLL; zw_printf('stage is renewed to %d\n', node_struct(ii).cw_stage); end
                                    node_struct(ii).backoffleft = int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
%                                     disp(CW_W(1)*2^(node_struct(ii).cw_stage));
                                    node_struct(ii).succeeded = 0;
                                    node_struct(ii).waitingchannel = 1;
                                else
                                    
                                    if node_struct(ii).inatim == 0 
                                        node_struct(ii).queue = node_struct(ii).queue - 1;
                                        node_struct(ii).success = node_struct(ii).success + 1;
                                        node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                    end
                                    temp = nearest_time - last_position;
                                    if node_struct(ii).inatim == 1 
                                        switch Actnode_cur(last_offset) 
                                            case  N_W/2
                                                packet_mactime(Round, 1) = temp*10;
                                            case  N_W/2-1
                                                packet_mactime(Round, 2) = temp*10;
                                            case  N_W/2-2
                                                packet_mactime(Round, 3) = temp*10;
                                            case  N_W/2-3
                                                packet_mactime(Round, 4) = temp*10;
                                            case  N_W/2-4
                                                packet_mactime(Round, 5) = temp*10;
                                            case  N_W/2-5
                                                packet_mactime(Round, 6) = temp*10;
                                            case  N_W/2-6
                                                packet_mactime(Round, 7) = temp*10;
                                            case  N_W/2-7
                                                packet_mactime(Round, 8) = temp*10;
                                            case  N_W/2-8
                                                packet_mactime(Round, 9) = temp*10;
                                        end
                                    end
                                    Actnode_cur(last_offset:last_offset+temp-1)=Actnode_cur(last_offset);
%                                     zw_printf('Actnode_cur(last_offset): %d\n', Actnode_cur(last_offset));
                                    Actnode_cur(last_offset+temp)=Actnode_cur(last_offset)-1;
                                    last_offset = last_offset + temp;
                                    last_position = nearest_time;
                                    node_struct(ii).next_event_time = Inf;
                                    node_struct(ii).next_event = EVENT_NONE;
                                    if DEBUG_COLL; zw_printf('[PSM:%d @ %d us] sending succeed\n', ii, nearest_time*10); end
                                    node_struct(ii).succeeded = 1;
                                    node_struct(ii).cw_stage = 0;
                                    
                                    
                                end
                               
                           end
                            
                    end

%                 end

            case WAKE_UP % wakeup event
%                 zw_printf('WAKEUP event\n');
                if node_struct(ii).type == DEVICE_LPL
                
                    if node_struct(ii).device == TRANSMITTER % transmitter
                        error(sprintf('[LPL:%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10));
    %                     if node_struct(ii).next_event == EVENT_SENT % is currently sending
    %                         node_struct(ii).sleep_event_time = node_struct(ii).next_event_time + 1;%ms to 10us
    %                         zw_printf('[%d @ %d] node is sending, defer to sleep at %d\n', ii, nearest_time, node_struct(ii).sleep_event_time);
    %                         
    %                     else
    %                         node_struct(ii).sleep_event_time = Inf;%ms to 10us
    %                         node_struct(ii).next_event = EVENT_NONE;
    %                         node_struct(ii).next_event_time = Inf;
    %                         zw_printf('[%d @ %d] node is not sending, sleep immediately\n', ii, nearest_time);
    %                         node_struct(ii).awake = 0;
    %                     end

                    else % receiver
                        if node_struct(ii).awake == 0
                            node_struct(ii).wakeup_event_time = node_struct(ii).wakeup_event_time + T_B*msto10us;%ms to 10us
                            if node_struct(ii).sleep_event_time == Inf; node_struct(ii).sleep_event_time = nearest_time; end
                            node_struct(ii).sleep_event_time = node_struct(ii).sleep_event_time + Tl_B*100;%ms to 10us
                            if DEBUG_WAKEUP; zw_printf('[LPL:%d @ %d us] wakes up, will go to sleep at %d us\n', ii, nearest_time*10, node_struct(ii).sleep_event_time*10); end
                            node_struct(ii).awake = 1;
                        else


                            if DEBUG_WAKEUP; error(sprintf('[LPL:%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10)); end
    %                         node_struct(ii).ducy_event_time = node_struct(ii).ducy_event_time + T*10;%*1e2;%ms to 10us
    %                         
    %                         zw_printf('[%d @ %d] goes to sleep, will wake up at %d\n', ii, nearest_time, node_struct(ii).ducy_event_time);
    %                         node_struct(ii).awake = 0;

                        end
                    end
                else %PSM, this is for ATIM window starting
                    node_struct(ii).succeeded = 0;
                    node_struct(ii).cw_stage = 0;
                    node_struct(ii).waitingchannel = 0;
                    node_struct(ii).inatim = 1;
                    node_struct(ii).backoffleft = -1;
                    node_struct(ii).last_awake = nearest_time;
                   
                    node_struct(ii).wakeup_event_time = Inf;%node_struct(ii).wakeup_event_time + int32(T_W*msto10us);%ms to 10us
                    if DEBUG_WAKEUP; zw_printf('[PSM:%d @ %d us] ATIM starts, ', ii, nearest_time*10); end
                    if node_struct(ii).sleep_event_time == Inf; node_struct(ii).sleep_event_time = nearest_time; end
                    node_struct(ii).sleep_event_time = node_struct(ii).sleep_event_time + int32(T_W*Rho_W*msto10us);%ms to 10us
                    if DEBUG_WAKEUP; zw_printf('will reach DATA at %d us, ', node_struct(ii).sleep_event_time*10); end
                    
%                     assert(node_struct(ii).queue >= 0);
                    
                    if node_struct(ii).queue > 0 || node_struct(ii).sat == 1%if having data in the queue
                        
                        if channel.idle == 1 %channel is idle
                            if node_struct(ii).next_event_time == Inf; node_struct(ii).next_event_time = nearest_time; end
                            node_struct(ii).next_event_time = node_struct(ii).next_event_time + DIFS;
                            node_struct(ii).backoffleft = int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
                            if DEBUG_WAKEUP; zw_printf('DIFS starts, backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
                            node_struct(ii).next_event = EVENT_BACKOFF;
%                             node_struct(ii).difsing = 1;
                        else
                            node_struct(ii).next_event = EVENT_NONE;
                            node_struct(ii).next_event_time = Inf;
                            if DEBUG_WAKEUP; zw_printf('channel is busy, waiting...\n'); end
                            node_struct(ii).waitingchannel = 1;
%                             if node_struct(ii).next_event_time == Inf; node_struct(ii).next_event_time = nearest_time; end
%                             node_struct(ii).next_event_time = node_struct(ii).next_event_time + channel.busyuntil - nearest_time + DIFS;
%                             node_struct(ii).backoffleft = int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
%                             if DEBUG_WAKEUP; zw_printf('backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
%                             node_struct(ii).next_event = EVENT_BACKOFF;
%                             node_struct(ii).difsing = 1;
                        end
                    else
                        node_struct(ii).sleepindw = 1;
                        if DEBUG_WAKEUP; zw_printf('no data to transmit\n'); end
                    end
                    
%                     node_struct(ii).waiting = 0;
%                     node_struct(ii).awake = 1;
%                     node_struct(ii).last_awake = nearest_time;
%                     node_struct(ii).queue = node_struct(ii).queue + 1;
%                     node_struct(ii).txstart = nearest_time;
                    

                    channel.atimstart = nearest_time;
                   
                
                end
                
            case GOTO_SLEEP % sleep event
%                 zw_printf('GOTO SLEEP event\n');
                if node_struct(ii).type == DEVICE_LPL

                    if node_struct(ii).device == TRANSMITTER % transmitter

                        if node_struct(ii).next_event == EVENT_SENT % is currently sending
                            node_struct(ii).sleep_event_time = node_struct(ii).next_event_time + 1;%ms to 10us
                            if DEBUG_SLEEP; zw_printf('[LPL:%d @ %d us] node is sending, defer to sleep at %d us\n', ii, nearest_time*10, node_struct(ii).sleep_event_time*10); end

                        else
                            node_struct(ii).queue = node_struct(ii).queue - 1;
                            node_struct(ii).count = node_struct(ii).count + 1;
                            node_struct(ii).drop = node_struct(ii).drop + 1;

%                             assert(node_struct(ii).queue >= 0);
                            if node_struct(ii).queue > 0 || node_struct(ii).sat == 1 %has a packet in the queue
                                node_struct(ii).sleep_event_time = nearest_time +1+ T_B*msto10us;%ms to 10us
                                node_struct(ii).next_event_time = nearest_time + 1;
                                node_struct(ii).next_event = EVENT_BACKOFF;
                                node_struct(ii).waiting = 0;
                                node_struct(ii).awake = 1;
                                node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                node_struct(ii).last_awake = nearest_time;
                                node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                                node_struct(ii).txstart = nearest_time;
                                if DEBUG_SLEEP; zw_printf('[LPL:%d @ %d us] node is not sending, and there is a new packet, drop and begin to backoff at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                            else
                                node_struct(ii).sleep_event_time = Inf;%ms to 10us
                                node_struct(ii).next_event = EVENT_NONE;
                                node_struct(ii).next_event_time = Inf;
                                if DEBUG_SLEEP; zw_printf('[LPL:%d @ %d us] node is not sending, and there is no new packet, drop and sleep immediately\n', ii, nearest_time*10); end
                                node_struct(ii).awake = 0;
                                node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                                node_struct(ii).txtime(node_struct(ii).count) = nearest_time - node_struct(ii).txstart;
                                node_struct(ii).waiting = 1;
                            end
                        end

                    else % receiver
                        if node_struct(ii).awake == 0
                            error(sprintf('[LPL:%d @ %d us] something is wrong if reaching here\n', ii, nearest_time*10));
    %                         node_struct(ii).ducy_event_time = node_struct(ii).ducy_event_time + T*10;%*1e2;%ms to 10us
    %                         zw_printf('[%d @ %d] wakes up, will go to sleep at %d\n', ii, nearest_time, node_struct(ii).ducy_event_time);
    %                         node_struct(ii).awake = 1;
                        else
                            node_struct(ii).sleep_event_time = Inf;

    %                         if DEBUG_SLEEP; zw_printf('[%d @ %d us] goes to sleep, will wake up at %d us\n', ii, nearest_time*10, node_struct(ii).wakeup_event_time*10); end
                            node_struct(ii).awake = 0;

                        end
                    end
                else %PSM, this is for DATA window starting
                    node_struct(ii).inatim = 0;
                    node_struct(ii).cw_stage = 0;
                    node_struct(ii).backoffleft = -1;
                    node_struct(ii).waitingchannel = 0;
                    node_struct(ii).sleep_event_time = Inf;
                    if node_struct(ii).wakeup_event_time == Inf; node_struct(ii).wakeup_event_time = nearest_time; end
                    node_struct(ii).wakeup_event_time = node_struct(ii).wakeup_event_time + int32(T_W*(1-Rho_W)*msto10us);%ms to 10us
                    if DEBUG_SLEEP; zw_printf('[PSM:%d @ %d us] DATA starts, node will reach next ATIM at %d us, ', ii, nearest_time*10, node_struct(ii).wakeup_event_time*10); end

                    if node_struct(ii).sleepindw == 1 %
                        
                        node_struct(ii).sleepindw = 0;
                        node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;
                        if DEBUG_SLEEP; zw_printf('no data, goes to sleep\n'); end
                       
                    else %check if transmission in atim succeed.
                        
                        if node_struct(ii).succeeded %atim packet is transmitted successfully
                            if DEBUG_SLEEP; zw_printf('ATIM succeed, will continue in DATA, '); end
                            if channel.idle == 1 %channel is idle
                            
                                if node_struct(ii).next_event_time == Inf; node_struct(ii).next_event_time = nearest_time; end
                                node_struct(ii).next_event_time = node_struct(ii).next_event_time + DIFS;
                                node_struct(ii).backoffleft = int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
                                if DEBUG_SLEEP; zw_printf('DIFS starts, backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
                                node_struct(ii).next_event = EVENT_BACKOFF;
%                                 node_struct(ii).difsing = 1;
                            else
                                node_struct(ii).next_event = EVENT_NONE;
                                node_struct(ii).next_event_time = Inf;
                                if DEBUG_WAKEUP; zw_printf('channel is busy, waiting...\n'); end
                                node_struct(ii).waitingchannel = 1;
                            end
                            
                        else
                        
                            if DEBUG_SLEEP; zw_printf('ATIM fails, goes to sleep\n'); end
                            node_struct(ii).next_event = EVENT_NONE;
                            node_struct(ii).next_event_time = Inf;
                            node_struct(ii).awaketime = node_struct(ii).awaketime + nearest_time - node_struct(ii).last_awake;

                        
                        end
                    end
                    
                    channel.datastart = nearest_time;

                    
                    
                    
                end
                
                

            case DATA_ARRIVAL %data arrival event
%                 zw_printf('DATA_ARRIVAL event\n');
                if node_struct(ii).type == DEVICE_LPL

                    if DEBUG_DR; zw_printf('[LPL:%d @ %d us] packet arrives..., ', ii, nearest_time*10); end
                    if node_struct(ii).waiting == 1
                        node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/TR_B)*log(myrand)*1e5);
                        if SEND_ONE_PACKET || node_struct(ii).sat == 1; node_struct(ii).dtar_event_time = Inf; end
                        if DEBUG_DR; zw_printf('next packet will arrive at %d us, ', node_struct(ii).dtar_event_time*10); end
                        if node_struct(ii).sleep_event_time == Inf; node_struct(ii).sleep_event_time = nearest_time; end
                        node_struct(ii).sleep_event_time = node_struct(ii).sleep_event_time +1 + T_B*msto10us;%ms to 10us
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
                        node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/TR_B)*log(myrand)*1e5);
                        if SEND_ONE_PACKET || node_struct(ii).sat == 1; node_struct(ii).dtar_event_time = Inf; end
                        if DEBUG_DR; zw_printf('next packet will arrive at %d us, ', node_struct(ii).dtar_event_time*10); end
                        if DEBUG_DR; zw_printf('node is sending, put the packet into the queue\n'); end
                        node_struct(ii).queue = node_struct(ii).queue + 1;
                    end
                else %PSM, simply put the data into the queue
                    if DEBUG_DR; zw_printf('[PSM:%d @ %d us] packet arrives..., ', ii, nearest_time*10); end
                    node_struct(ii).dtar_event_time = node_struct(ii).dtar_event_time + ceil((-1/TR_W)*log(myrand)*1e5);
                    if SEND_ONE_PACKET; node_struct(ii).dtar_event_time = Inf; end
                    if DEBUG_DR; zw_printf('next packet will arrive at %d us, ', node_struct(ii).dtar_event_time*10); end
                    if DEBUG_DR; zw_printf('the packet is put in the queue\n'); end
                    node_struct(ii).queue = node_struct(ii).queue + 1;
                end

        end    
    
    
    end
    
    transmittings = find([node_struct.next_event] == EVENT_SENT);
%     disp(transmittings);
    if isempty(transmittings)
        if channel.idle == 0
            channel.idle = 1;
            if DEBUG_CHANNEL; zw_printf('[Channel @ %d us] channel becomes idle\n', nearest_time*10); end
            
            for ii=1:N_W/2 %only considering the transmitters for PSM

                if node_struct(ii).waitingchannel == 1
    %                 node_struct(ii).next_event_time = node_struct(ii).next_event_time + channel.busyfor + DIFS;
    %                 if 1; zw_printf('backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
    %                 node_struct(ii).next_event = EVENT_BACKOFF;                
                    assert(node_struct(ii).next_event_time == Inf);
                    if node_struct(ii).next_event_time == Inf; node_struct(ii).next_event_time = nearest_time; end
                    node_struct(ii).next_event_time = node_struct(ii).next_event_time + DIFS;
                    if node_struct(ii).backoffleft == -1
                        node_struct(ii).backoffleft = int32(myrand*CW_W(1)*2^(node_struct(ii).cw_stage));
                    end
                    if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] DIFS starts, backoff will begin at %d us\n', ii, nearest_time*10, node_struct(ii).next_event_time*10); end
                    node_struct(ii).next_event = EVENT_BACKOFF;
                   
                    
                    node_struct(ii).waitingchannel = 0;
                    

%                 elseif node_struct(ii).backingoff == 1
%                     node_struct(ii).next_event = EVENT_NONE;
%                     node_struct(ii).next_event_time = Inf;
% %                     node_struct(ii).backoffleft = node_struct(ii).next_event_time - nearest_time;
% %                     node_struct(ii).next_event_time = nearest_time + channel.busyfor + DIFS;
% %                     if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] freezing for %d us and backoff will restart at %d us\n', ii, nearest_time*10, channel.busyfor*10, node_struct(ii).next_event_time*10); end
% %                     node_struct(ii).next_event = EVENT_BACKOFF;
% %                     node_struct(ii).difsing = 1;
% %                     node_struct(ii).backingoff = 0;

                end

            end               
            
        end
    else
        if channel.idle == 1
            channel.idle = 0;
            if DEBUG_CHANNEL; zw_printf('[Channel @ %d us] channel becomes busy...', nearest_time*10); end
            if size(transmittings,2) > 1 %multiple transmitings.
                if DEBUG_CHANNEL; zw_printf('and there is a collision...'); end
                channel.collide = 1;
                if any(transmittings<=N_W/2) && any(transmittings>N_W/2)
                    if DEBUG_CHANNEL; zw_printf('which is a inter-collision\n'); end
%                     disp(nearest_time*10);
%                     disp(transmittings);
                else
                    if DEBUG_CHANNEL; zw_printf('\n'); end
                end
            else
                if DEBUG_CHANNEL; zw_printf('and transmission is collision free\n'); end
                channel.collide = 0;
            end
            
            for ii=1:N_W/2 %only considering the transmitters for PSM

                if node_struct(ii).next_event == EVENT_BACKOFF
    %                 node_struct(ii).next_event_time = node_struct(ii).next_event_time + channel.busyfor + DIFS;
    %                 if 1; zw_printf('backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
    %                 node_struct(ii).next_event = EVENT_BACKOFF;                
                    node_struct(ii).next_event = EVENT_NONE;
                    node_struct(ii).next_event_time = Inf;
                    node_struct(ii).waitingchannel = 1;
                    if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] freezing until window ends or channel idle\n', ii, nearest_time*10); end


                elseif node_struct(ii).next_event == EVENT_SENSE
                    node_struct(ii).backoffleft = node_struct(ii).next_event_time - nearest_time;
                    node_struct(ii).next_event = EVENT_NONE;
                    node_struct(ii).next_event_time = Inf;
                    node_struct(ii).waitingchannel = 1;
                    if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] freezing until window ends or channel idle\n', ii, nearest_time*10); end
%                     node_struct(ii).next_event_time = nearest_time + channel.busyfor + DIFS;
%                     if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] freezing for %d us and backoff will restart at %d us\n', ii, nearest_time*10, channel.busyfor*10, node_struct(ii).next_event_time*10); end
%                     node_struct(ii).next_event = EVENT_BACKOFF;
%                     node_struct(ii).difsing = 1;
%                     node_struct(ii).backingoff = 0;

                end

            end            
            
            
        end
        
    end
    
    if channel.atimstart == nearest_time
        
        if Round > 0
            if floor(Round/print_density) > current_round;
              zw_printf('[%d us] Round: %d\n', nearest_time*10, Round);
              current_round = floor(Round/print_density);
            end
            
            temp = nearest_time - last_position;
            Actnode_cur(last_offset:last_offset+temp)=Actnode_cur(last_offset);
            Actnode_stat(1:T_W*1e3/timeslot) = Actnode_stat(1:T_W*1e3/timeslot) + Actnode_cur(1:T_W*1e3/timeslot);
            
            
%             last_offset = last_offset + temp;
%             last_position = nearest_time;
            if showdetailed_nt
                if singlefig
                    figure(afig);
                else
                    figure;
                end

                plot(Actnode_cur, 'DisplayName','Actnode_cur', 'Color',col(Round,:));
                xlabel('Time');ylabel('Average # Active nodes A(t)');
                title(sprintf('N=%d, W_0=%d, T_W=%d, d=%d, Drawing (%d/%d) round', N_W/2, cw_W, T_W, L_W, Round, sim_time/T_W));
                hold on;
            end
        end
        
        Actnode_cur(1:T_W*1e3/timeslot)=-1;
        Actnode_cur(1)=N_W/2;
        last_offset = 1;
        last_position = nearest_time;

        Round = Round + 1;        
        
    end
    
    if channel.datastart == nearest_time
        temp = nearest_time - last_position;
        Actnode_cur(last_offset:last_offset+temp-1)=Actnode_cur(last_offset);
        Actnode_cur(last_offset+temp)=N_W/2-Actnode_cur(last_offset);
        packet_sucinatim(:,Round) = Actnode_cur(last_offset+temp);
        last_offset = last_offset + temp;
        last_position = nearest_time;
        
    
    end
    
%     if channel.idle2busy==1% && channel.collide == 1 %if channel becomes busy at this time instance, then we need to tell PSM nodes about this
%         if DEBUG_CHANNEL; zw_printf('[Channel @ %d us] channel busy for %d us\n', nearest_time*10, channel.busyfor*10); end
%         
%         for ii=1:N_W/2 %only considering the transmitters for PSM
%         
%             if node_struct(ii).difsing == 1
% %                 node_struct(ii).next_event_time = node_struct(ii).next_event_time + channel.busyfor + DIFS;
% %                 if 1; zw_printf('backoff will begin at %d us\n', node_struct(ii).next_event_time*10); end
% %                 node_struct(ii).next_event = EVENT_BACKOFF;                
%                 
%             elseif node_struct(ii).backingoff == 1
%                 node_struct(ii).backoffleft = node_struct(ii).next_event_time - nearest_time;
%                 node_struct(ii).next_event_time = nearest_time + channel.busyfor + DIFS;
%                 if DEBUG_CHANNEL; zw_printf('[PSM:%d @ %d us] freezing for %d us and backoff will restart at %d us\n', ii, nearest_time*10, channel.busyfor*10, node_struct(ii).next_event_time*10); end
%                 node_struct(ii).next_event = EVENT_BACKOFF;
%                 node_struct(ii).difsing = 1;
%                 node_struct(ii).backingoff = 0;
%                 
%             end
%             
%         end
%         
%         channel.idle2busy = 0;
%         channel.busyfor = 0;
%     end
    
    
    
    
    t = nearest_time;
    
    if t>=true_time
        
       break;
       
    end
    
    
end



%             last_offset = last_offset + temp;
%             last_position = nearest_time;

asfig = figure;
plot(Actnode_stat/(Round-1), 'DisplayName','Actnode_cur');
xlabel('Time');ylabel('Average # Active nodes A(t)');
% title(sprintf('N=%d, W_0=%d, T_W=%d, d=%d', N_W/2, cw_W, T_W, L_W));
title(sprintf('N_W=%d, W_W=%d, T_W=%d us, L_W=%d, N_B=%d, W_B=%d, T_B=%d us, L_B=%d', ...
                N_W/2, cw_W, T_W, L_W, N_B/2, cw_B, T_B, L_B));



toc;
if DEBUG_LOG
    diary off;
end

% outname=sprintf('stat(%d-%d-%d-%g-%d-%g).mat', sim_time/1000, N_B, T_B, Rho_B, CW_B(1), TR_B);
% 
% save(outname, 'node_struct');

out = node_struct;
% out(2) = packet_mactime*timeslot;

end
