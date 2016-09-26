#!/bin/bash

Num_files=$(ls log*.out|wc -l)

for j in $(seq 1 $Num_files)
do
log_file=$(ls .|awk "NR==$j")
Wifi_Nodes=$(grep 'WIFI_NODES:' $log_file| awk -F: '{print $2}')
Wifi_Nodes=$(echo "$Wifi_Nodes / 2"|bc)
Avg_Coll_rate_per_node=0
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
done

Wifi_Pckts=$(grep 'Success to transmit' $log_file|wc -l)
Wifi_Pckt_size=$(grep 'WIFI_PCKT_SIZE:' $log_file| awk -F: '{print $2}')
Zigbee_Nodes=$(grep 'ZIGBEE_NODES:' $log_file| awk -F: '{print $2}')
Zigbee_Pckts=$(grep 'Received' $log_file|wc -l)
Zigbee_Pckt_size=$(grep 'ZIGBEE_PCKT_SIZE:' $log_file| awk -F: '{print $2}')
Time_duration=$(grep 'DURATION:' $log_file| awk -F: '{print $2}')

Wifi_Throughput=$(echo "scale=3; $Wifi_Pckt_size * $Wifi_Pckts * 8 / $Time_duration / 54000000"|bc)
Zigbee_Throughput=$(echo "scale=3; $Zigbee_Pckt_size * $Zigbee_Pckts * 8 / $Time_duration / 250000"|bc)
Avg_Coll_rate_per_node=$(echo "scale=3; $Avg_Coll_rate_per_node/ $Wifi_Nodes"|bc)

echo "\n-----PARAMETERS ------\n "
grep 'WIFI' $log_file
grep 'ZIGBEE' $log_file
echo "\nNo. of Wifi Packets Successfully Transmitted: " $Wifi_Pckts
echo "No. of Zigbee Packets Successfully Transmitted: " $Zigbee_Pckts
echo "Wifi Throughput: " $Wifi_Throughput "bps"
echo "Zigbee Throughput: " $Zigbee_Throughput "bps" 
echo "Average Collision rate per node: " $Avg_Coll_rate_per_node 
done


