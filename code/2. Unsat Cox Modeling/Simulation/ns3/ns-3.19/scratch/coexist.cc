/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014 Texas A&M University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as 
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Author: Raghavan S V <raghavan@tamu.edu> 
 */


#include "ns3/log.h"
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/lr-wpan-module.h"
#include "ns3/propagation-loss-model.h"
#include "ns3/propagation-delay-model.h"
#include "ns3/simulator.h"
#include <ns3/multi-model-spectrum-channel.h>
#include <ns3/constant-position-mobility-model.h>
#include "ns3/sw-module.h"
#include "ns3/sw-mac-csma-helper.h"
#include "ns3/sw-phy-basic-helper.h"
#include <vector>
#include <iostream>

using namespace ns3;

int LrWpan_Nodes = 6;
int Wifi_Nodes = 6;
int numPackets = 10000;
int Wifi_pkt = 1600;
int LrWpan_pkt = 88;
char Addr[256][6];
int Wifi_arrival_time = 600;
int LrWpan_arrival_time = 1000;

NS_LOG_COMPONENT_DEFINE("Coexistence");

void PopulateAddrs ()
{
  int i,j;
  j=1;
  Addr[0][0] = '0';
  Addr[0][1] = '0';
  Addr[0][2] = ':';
  Addr[0][3] = '0';
  Addr[0][4] = '0';
  Addr[0][5] = '\0';
  while(j<256)
  {
        for(i=0;i<4;i++)
        {
                Addr[j][i] = Addr[j-1][i];
        }    
        if((Addr[j-1][i] <= '8') || ((Addr[j-1][i] >= 'a') && (Addr[j-1][i] <= 'e')))
        {
                Addr[j][i] = Addr[j-1][i]+1;
        }
        else if(Addr[j-1][i] == '9')
        {
                Addr[j][i] = 'a';
        }
        else
        {
                Addr[j][i] = '0';
                if((Addr[j-1][i-1] <= '8') || ((Addr[j-1][i-1] >= 'a') && (Addr[j-1][i-1] <= 'e')))
                {
                        Addr[j][i-1]++;
                }
                else if(Addr[j-1][i-1] == '9')
                {
                        Addr[j][i-1] = 'a';
                }
                else
                {
                        Addr[j][i-1] = '0';     
                }                                   
        }
        Addr[j][5] = '\0';
        j++;
  }
}

#if 1
void GenerateLrWpanTraffic (NetDeviceContainer devices, int pck_count, McpsDataRequestParams params, int j)
{

  Ptr<Packet> p = Create<Packet> (LrWpan_pkt);
  params.m_dstAddr = Mac16Address (Addr[j]);
  Ptr<LrWpanMac> bmac = devices.Get(j-LrWpan_Nodes/2)->GetObject<LrWpanNetDevice> ()->GetMac ();
  bmac->McpsDataRequest(params, p);
  Simulator::Schedule (MicroSeconds (LrWpan_arrival_time), &GenerateLrWpanTraffic, devices, pck_count-1, params, j);
}

void GenerateWifiTraffic (NetDeviceContainer devices, int pck_count, int k)
{

  Ptr<Packet> p = Create<Packet> (Wifi_pkt);  
  devices.Get(k - Wifi_Nodes/2)->Send(p, devices.Get(k)->GetAddress (), 1);
  Simulator::Schedule (MicroSeconds (Wifi_arrival_time), &GenerateWifiTraffic, devices, pck_count-1, k);
}
#else
void GenerateWifiTraffic (NetDeviceContainer devices, Ptr<Packet> p, int k)
{
  devices.Get(k - Wifi_Nodes/2)->Send(p, devices.Get(k)->GetAddress (), 1);
}
#endif 
static void DataIndication (McpsDataIndicationParams params, Ptr<Packet> p, uint32_t m_pktNumSuccess)
{
  NS_LOG_UNCOND (Simulator::Now()<< " Node: "<<params.m_dstAddr<<" Received Packet "<< m_pktNumSuccess <<" from Node "<< params.m_srcAddr);
}

int main (int argc, char *argv[])
{
  /* ------ 802.15.4 ----------- */
 
  bool verbose = false;
  int Stop_Time = 2;
  int CW_tune_LrWpan = 80;
  int CW_tune_Wifi = 16;
  CommandLine cmd;

  cmd.AddValue ("verbose", "turn on all log components", verbose);
  cmd.AddValue ("LrWpan_Nodes", "Number of 802.15.4 nodes", LrWpan_Nodes);
  cmd.AddValue ("Wifi_Nodes", "Number of Wifi nodes", Wifi_Nodes);
  cmd.AddValue ("numPackets", "Number of packets", numPackets);
  cmd.AddValue ("Stop_Time", "Stop time for simulation", Stop_Time);
  cmd.AddValue ("CW_tune_LrWpan", "Contention Window Size tuning - Zigbee", CW_tune_LrWpan);
  cmd.AddValue ("CW_tune_Wifi", "Contention Window Size tuning - Wifi", CW_tune_Wifi);
  cmd.AddValue ("LrWpan_arrival_time", "Interarrival Time - LrWpan", LrWpan_arrival_time);
  cmd.AddValue ("Wifi_arrival_time", "Interarrival Time - Wifi", Wifi_arrival_time);


  cmd.Parse (argc, argv);

  NS_LOG_UNCOND ("---------- PARAMETERS ------------");
  NS_LOG_UNCOND ("WIFI_NODES:"<<Wifi_Nodes<<"\nZIGBEE_NODES:"<<LrWpan_Nodes<<"\nDURATION:"<<(Stop_Time-1)<<"\nWIFI_CW_SIZE:"<< CW_tune_Wifi<<"\nZIGBEE_CW_SIZE:"<<CW_tune_LrWpan<<"\nWIFI_PCKT_SIZE:"<<Wifi_pkt<<"\nZIGBEE_PCKT_SIZE:"
<<LrWpan_pkt<<"\nWIFI_ARRIVAL_TIME:"<<Wifi_arrival_time<<"\nZIGBEE_ARRIVAL_TIME:"<<LrWpan_arrival_time
<<"\n");

  LrWpanHelper lrWpanHelper;
  if (verbose)
    {
        LogComponentEnable ("LrWpanMac", LOG_LEVEL_DEBUG);
        LogComponentEnable ("LrWpanPhy", LOG_LEVEL_DEBUG);
        LogComponentEnable ("LrWpanCsmaCa", LOG_LEVEL_DEBUG);
        LogComponentEnable ("SwMacCsma", LOG_LEVEL_DEBUG);
        LogComponentEnable ("SwPhy", LOG_LEVEL_DEBUG);
    }
  
  
  PopulateAddrs ();
  NodeContainer lrwpan_nodes;
  lrwpan_nodes.Create (LrWpan_Nodes);

  LrWpanHelper lrwpan;
  NetDeviceContainer lrwpan_devices = lrwpan.Install (lrwpan_nodes);
  
  MobilityHelper lrwpan_mobility;
  Ptr<ListPositionAllocator> positionAlloc_lrwpan = CreateObject<ListPositionAllocator> ();

  McpsDataIndicationCallback Icb;
  Icb = MakeCallback (&DataIndication);

  for(int i=0; i<LrWpan_Nodes; i++)
  {
        lrwpan_devices.Get(i)->GetObject<LrWpanNetDevice> ()->SetAddress (Mac16Address (Addr[i]));
        positionAlloc_lrwpan->Add (Vector (0.0, 0.0, 0.0));
        lrwpan_devices.Get(i)->GetObject<LrWpanNetDevice> ()->GetMac ()->SetMcpsDataIndicationCallback (Icb); 
	lrwpan_devices.Get(i)->GetObject<LrWpanNetDevice> ()->GetMac ()->OS_Delay_lrwpan = MicroSeconds(6600);   
	lrwpan_devices.Get(i)->GetObject<LrWpanNetDevice> ()->GetCsmaCa ()->setMacMaxBE (CW_tune_LrWpan);   
  }

  lrwpan_mobility.SetPositionAllocator (positionAlloc_lrwpan);
  lrwpan_mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  lrwpan_mobility.Install (lrwpan_nodes);

  McpsDataRequestParams params;
  params.m_srcAddrMode = 2;
  params.m_dstAddrMode = 2;
  params.m_dstPanId = 0;
  params.m_msduHandle = 0;
  params.m_txOptions = 0;

#if 1

  for(int j=LrWpan_Nodes/2;j<LrWpan_Nodes;j++)
  {
        Simulator::Schedule(Seconds(1), &GenerateLrWpanTraffic, lrwpan_devices, numPackets, params, j);    
  }

#else
  
  for(int i=0;i<numPackets;i++)
  {
        for(int j=LrWpan_Nodes/2;j<LrWpan_Nodes;j++)
        {
          Ptr<Packet> p = Create<Packet> (LrWpan_pkt);
          params.m_dstAddr = Mac16Address (Addr[j]);
          Simulator::Schedule (MilliSeconds (1000+i*10), &LrWpanMac::McpsDataRequest, lrwpan_devices.Get(j-LrWpan_Nodes/2)->GetObject<LrWpanNetDevice> ()->GetMac (), params, p);
        }
  }
#endif

  /*--------802.11 ----------- */
  
  NodeContainer nodes;
  nodes.Create (Wifi_Nodes);
  
  SwMacCsmaHelper swMac = SwMacCsmaHelper::Default ();
  SwPhyBasicHelper swPhy = SwPhyBasicHelper::Default ();
  SwHelper sw;
  NetDeviceContainer wifi_devices = sw.Install (nodes, lrwpan_devices.Get(0)->GetObject<LrWpanNetDevice> ()->DoGetChannel(), swPhy, swMac);

  MobilityHelper mobility;
  Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
  for(int i=0; i<Wifi_Nodes; i++)
  {
  positionAlloc->Add (Vector (0.0, 0.0, 0.0));
  wifi_devices.Get(i)->GetObject<SwNetDevice> ()-> GetMac () ->SetCwMin (CW_tune_Wifi);
  }
  mobility.SetPositionAllocator (positionAlloc);
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobility.Install (nodes);

#if 1

  for(int k=Wifi_Nodes/2; k<Wifi_Nodes; k++)    
  {
       Simulator::Schedule(Seconds(1),&GenerateWifiTraffic, wifi_devices, numPackets, k); 
  } 
    
#else   
  for(int i=0;i<numPackets;i++)
  {
        for(int k=Wifi_Nodes/2; k<Wifi_Nodes; k++)    
        {
        Ptr<Packet> p = Create<Packet> (Wifi_pkt);  
        Simulator::Schedule (MilliSeconds (1000+i), &GenerateWifiTraffic, wifi_devices, p, k);
        }
  }
#endif     
        
  Simulator::Stop (Seconds(Stop_Time));
  Simulator::Run ();
  Simulator::Destroy ();

  return 0;
}

