#!/bin/bash

Num_files=$(ls logs_sim/"$1"/log*.out|wc -l)

for j in $(seq 1 $Num_files)
do
log_file=$(ls logs_sim/"$1"/log*.out|awk "NR==$j")
Wifi_Nodes=$(grep 'WIFI_NODES:' $log_file| awk -F: '{print $2}')
Wifi_Nodes=$(echo "$Wifi_Nodes / 2"|bc)
Zigbee_Nodes=$(grep 'ZIGBEE_NODES:' $log_file| awk -F: '{print $2}')
Zigbee_Nodes=$(echo "$Zigbee_Nodes / 2"|bc)

Wifi_Pckts=$(grep 'Success to transmit' $log_file|wc -l)
Wifi_Pckt_size=$(grep 'WIFI_PCKT_SIZE:' $log_file| awk -F: '{print $2}')
Zigbee_Pckts=$(grep 'Received' $log_file|wc -l)
Zigbee_Pckt_size=$(grep 'ZIGBEE_PCKT_SIZE:' $log_file| awk -F: '{print $2}')
Time_duration=$(grep 'DURATION:' $log_file| awk -F: '{print $2}')

#Wifi stuff
Avg_Coll_rate_per_node=0
Wifi_avg_pkt_delay_per_node=0
Wifi_avg_queue_size_per_node=0
for i in $(seq 1 $Wifi_Nodes);
do
if [ $i -lt 16 ];
then
 i=$(echo "obase=16;$i"|bc)
 i=$(echo $i | tr '[:upper:]' '[:lower:]')
 i='0'$i
else
 i=$(echo "obase=16;$i"|bc)
 i=$(echo $i | tr '[:upper:]' '[:lower:]')
fi
Success=$(grep "Node Address : 00:00:00:00:00:$i" $log_file|wc -l)
Collision=$(grep "Packet transmission aborted - COLLISION (Wifi) Node: 00:00:00:00:00:$i" $log_file|wc -l)
Coll_rate_node=$(echo "scale=3; $Collision / ($Success + $Collision)"|bc)
Avg_Coll_rate_per_node=$(echo "scale=3; $Avg_Coll_rate_per_node + $Coll_rate_node"|bc)
Cumm_pkt_delay=$(grep "at Wifi Node 00:00:00:00:00:$i" $log_file| tail -1| awk '{print $14}'| sed 's/[a-zA-Z+]//g')
Cumm_queue_size=$(grep "at Wifi Node 00:00:00:00:00:$i" $log_file| tail -1| awk '{print $18}'| sed 's/[a-zA-Z+]//g')
Queue_chks=$(grep "at Wifi Node 00:00:00:00:00:$i" $log_file|wc -l);
Wifi_avg_pkt_delay_per_node=$(echo "scale=3; $Wifi_avg_pkt_delay_per_node + ($Cumm_pkt_delay/$Success)"|bc)
Wifi_avg_queue_size_per_node=$(echo "scale=3; $Wifi_avg_queue_size_per_node + ($Cumm_queue_size/$Queue_chks)"|bc)
done

#Zigbee stuff
Zigbee_avg_pkt_delay_per_node=0
Zigbee_avg_queue_size_per_node=0
for j in $(seq 1 $Zigbee_Nodes);
do
 j=`expr $j - 1`
if [ $j -lt 16 ];
then
 j=$(echo "obase=16;$j"|bc)
 j=$(echo $j | tr '[:upper:]' '[:lower:]')
 j='0'$j
else
 j=$(echo "obase=16;$j"|bc)
 j=$(echo $j | tr '[:upper:]' '[:lower:]')
fi
Success=$(grep "from Node 00:$j" $log_file|wc -l)
Cumm_pkt_delay=$(grep "at Zigbee Node 00:$j" $log_file| tail -1| awk '{print $14}'| sed 's/[a-zA-Z+]//g')
Cumm_queue_size=$(grep "at Zigbee Node 00:$j" $log_file| tail -1| awk '{print $18}'| sed 's/[a-zA-Z+]//g')
Queue_chks=$(grep "at Zigbee Node 00:$j" $log_file|wc -l);
Zigbee_avg_pkt_delay_per_node=$(echo "scale=3; $Zigbee_avg_pkt_delay_per_node + ($Cumm_pkt_delay/$Success)"|bc)
Zigbee_avg_queue_size_per_node=$(echo "scale=3; $Zigbee_avg_queue_size_per_node + ($Cumm_queue_size/$Queue_chks)"|bc)
done


Wifi_Throughput=$(echo "scale=3; $Wifi_Pckt_size * $Wifi_Pckts * 8 / $Time_duration / 54000000"|bc)
Wifi_avg_pckt_delay=$(echo "scale=3; $Wifi_avg_pkt_delay_per_node / $Wifi_Nodes"|bc)
Wifi_avg_queue_size=$(echo "scale=3; $Wifi_avg_queue_size_per_node / $Wifi_Nodes"|bc)
Zigbee_Throughput=$(echo "scale=3; $Zigbee_Pckt_size * $Zigbee_Pckts * 8 / $Time_duration / 250000"|bc)
Zigbee_avg_pckt_delay=$(echo "scale=3; $Zigbee_avg_pkt_delay_per_node / $Zigbee_Nodes"|bc)
Zigbee_avg_queue_size=$(echo "scale=3; $Zigbee_avg_queue_size_per_node / $Zigbee_Nodes"|bc)

Avg_Coll_rate_per_node=$(echo "scale=3; $Avg_Coll_rate_per_node / $Wifi_Nodes"|bc)

echo "-----PARAMETERS ------"
grep 'WIFI' $log_file
grep 'ZIGBEE' $log_file
echo "No. of Wifi Packets Successfully Transmitted: " $Wifi_Pckts
echo "No. of Zigbee Packets Successfully Transmitted: " $Zigbee_Pckts
echo "Wifi Throughput: " $Wifi_Throughput "bps"
echo "Zigbee Throughput: " $Zigbee_Throughput "bps" 
echo "Wifi Average MAC Service Packet Delay: " $Wifi_avg_pckt_delay "ns" 
echo "Zigbee Average MAC Service Packet Delay: " $Zigbee_avg_pckt_delay "ns" 
echo "Wifi Average Queue Size: " $Wifi_avg_queue_size
echo "Zigbee Average Queue Size: " $Zigbee_avg_queue_size  
echo "Wifi Average Collision rate per node: " $Avg_Coll_rate_per_node 
done


