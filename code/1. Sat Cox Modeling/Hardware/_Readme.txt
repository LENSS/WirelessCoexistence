Here are the code for the hardware based experiment of coexisting 4 TelosB motes and 2 Mikrotik routers.
The original code for TelosB is missing, but since it is very similar to the sample code in Tinyos called TxThroughput, I put it here.
Basically it just continuously transmits data to measure the maximum throughput, just like what we did for saturation scenario.
In terms of the WiFi side, we simply use the tool iperf to continuously send UDP packets.
