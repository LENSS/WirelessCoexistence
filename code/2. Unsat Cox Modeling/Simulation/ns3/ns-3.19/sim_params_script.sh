#!/bin/bash 

verbose=false
Wifi_pckt_size=1600
Zigbee_pckt_size=100
CW_min_Wifi=16
CW_max_Wifi=1024
CW_min_zigbee=320
CW_max_zigbee=80
Wifi_arrival_time=5000
LrWpan_arrival_time=1000
Stop_Time=3
Wifi_Nodes=10
Zigbee_Nodes=20
i=0
foldername="$(date +"%m-%d-[%T]")"
mkdir logs_sim/$foldername
for Wifi_pckt_size in 1600
do
 for CW_min_Wifi in 16 
 do
  for CW_max_Wifi in 1024
  do
   for Zigbee_pckt_size in 88
   do
    for CW_max_zigbee in 70 
    do
     for Wifi_Nodes in $Wifi_Nodes
     do
      for Zigbee_Nodes in $Zigbee_Nodes
      do
       for Wifi_arrival_time in $Wifi_arrival_time
       do
	for LrWpan_arrival_time in $LrWpan_arrival_time #(seq 200000 100000 600000) #middle value is the step-interval 
	do        
	 i=$(echo "$i + 1"|bc)
         filename="logs_sim/$foldername/log_$i-$Zigbee_Nodes-$Wifi_Nodes-$CW_max_zigbee-$CW_min_Wifi-$Wifi_arrival_time-$LrWpan_arrival_time.out"
         ./waf --run "scratch/coexist --verbose=$verbose --LrWpan_Nodes=$Zigbee_Nodes --Wifi_Nodes=$Wifi_Nodes --Stop_Time=$Stop_Time --CW_tune_LrWpan=$CW_max_zigbee --CW_tune_Wifi=$CW_min_Wifi --LrWpan_arrival_time=$LrWpan_arrival_time --Wifi_arrival_time=$Wifi_arrival_time" >$filename 2>&1 
	done
       done      
      done
     done
    done
   done
  done
 done
done

#Running script to generate summary
sh logs_sim/script_log.sh $foldername >logs_sim/$foldername/summary.txt



