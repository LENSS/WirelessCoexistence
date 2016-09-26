#!/bin/bash

Values=$(grep 'PARAMETERS' $1|wc -l)
filename="data_plot.txt"

echo '"Contention Window Sizes - WIFI" "ZIGBEE_CW_SIZE" "802.11-throughput-NS3" "BoX-MAC througput-NS3" "Analytical-throughput-802.11" "Analytical throughput-Bmac"'>$filename

for i in $(seq 1 $Values)
do
#grep 'WIFI_CW_SIZE:' $1 | awk "NR==$i" | awk -F: '{print $2}'  | tr '\n' '\t' >>$filename
echo $i | tr '\n' '\t' >>$filename
grep 'ZIGBEE_CW_SIZE:' $1 | awk "NR==$i" | awk -F: '{print $2}'  | tr '\n' '\t' >>$filename
grep 'Wifi Throughput: ' $1 | awk "NR==$i" | awk -F: '{print $2}'  | tr -d 'bps,' | tr '\n' '\t' >>$filename
grep 'Zigbee Throughput: ' $1 | awk "NR==$i" | awk -F: '{print $2}' | tr -d 'bps,' | tr '\r' '\t' >>$filename
#grep "aggrThr_trp_w =" $2| awk "NR==$i" | sed 's/[a-zA-Z_=\n]//g' | tr '\n\r' '\t' >> $filename
#grep "aggrThr_trp_b =" $2| awk "NR==$i" | sed 's/[a-zA-Z_=]//g' | tr '\r' '\t' >> $filename
done
