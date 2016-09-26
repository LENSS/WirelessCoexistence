/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2010 University of Arizona
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
 * Author: Junseok Kim <junseok@email.arizona.edu> <engr.arizona.edu/~junseok>
 */

#include "ns3/simulator.h"
#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "ns3/mac48-address.h"
#include "sw-mac.h"
#include "sw-phy.h"

NS_LOG_COMPONENT_DEFINE ("SwPhy");

namespace ns3 {

NS_OBJECT_ENSURE_REGISTERED (SwPhy);

TypeId
SwPhy::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SwPhy")
    .SetParent<Object> ()
    .AddConstructor<SwPhy> ()
    .AddAttribute ("PreambleDuration",
                   "Duration (us) of Preamble of PHY Layer",
                   TimeValue (MicroSeconds (16)),
                   MakeTimeAccessor (&SwPhy::m_preambleDuration),
                   MakeTimeChecker ())
    .AddAttribute ("TrailerSize",
                   "Size of Trailer (e.g. FCS) (bytes)",
                   UintegerValue (2),
                   MakeUintegerAccessor (&SwPhy::m_trailerSize),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("HeaderSize",
                   "Size of Header (bytes)",
                   UintegerValue (3),
                   MakeUintegerAccessor (&SwPhy::m_headerSize),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("SinrTh",
                   "SINR Threshold",
                   DoubleValue (2),
                   MakeDoubleAccessor (&SwPhy::m_sinrTh),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("CsPowerTh",
                   "Carrier Sense Threshold (dBm)",
                   DoubleValue (-110),
                   MakeDoubleAccessor (&SwPhy::m_csTh),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("TxPower",
                   "Transmission Power (dBm)",
                   DoubleValue (10),
                   MakeDoubleAccessor (&SwPhy::SetTxPower),
                   MakeDoubleChecker<double> ())
    ;
  return tid;
}

SwPhy::SwPhy ()
 : m_device (0),
   m_mac (0),
   m_channel (0),
   m_rxPsd(0),
   m_EndRxRequest (),
   m_EndCTRXRequest ()
{
  m_rxTotalPower = 0.0;
  m_pktRx = 0;
  SwSpectrumValueHelper psdHelper;
  m_txPsd = psdHelper.CreateTxPowerSpectralDensity (10, 11);
  m_noise = psdHelper.CreateNoisePowerSpectralDensity (11);
  m_csBusy = false;
  m_csBusyEnd = Seconds (0);
  txDuration = Seconds (0);
  duration = Seconds (0);
}
SwPhy::~SwPhy ()
{
  Clear ();
}
void
SwPhy::Clear ()
{
  m_pktRx = 0;
}

Ptr<MobilityModel>
SwPhy::GetMobility ()
{
  NS_LOG_FUNCTION (this);
  return m_mobility;
}

Ptr<NetDevice>
SwPhy::GetDevice ()
{
  NS_LOG_FUNCTION (this);
  return m_device;
}

void
SwPhy::SetDevice (Ptr<NetDevice> device)
{
  m_device = device;
}

void
SwPhy::SetMobility (Ptr<MobilityModel> m)
{
  NS_LOG_FUNCTION (this << m);
  m_mobility = m;
}

void
SwPhy::SetMac (Ptr<SwMac> mac)
{
  m_mac = mac;
}
void
SwPhy::SetChannel (Ptr<SpectrumChannel> channel)
{
  m_channel = channel;
}
void
SwPhy::SetTxPower (double dBm)
{
  m_txPower = dBm;
}

//-----------------------------------------------------------------
Ptr<SpectrumChannel>
SwPhy::GetChannel ()
{
  return m_channel;
}

Ptr<const SpectrumModel>
SwPhy::GetRxSpectrumModel () const
{
  if (m_txPsd)
    {
      return m_txPsd->GetSpectrumModel ();
    }
  else
    {
      return 0;
    }
}

Mac48Address
SwPhy::GetAddress ()
{
  return m_mac->GetAddress ();
}
double
SwPhy::GetTxPower ()
{
  return m_txPower;
}

Ptr<AntennaModel>
SwPhy::GetRxAntenna ()
{
  return m_antenna;
}

void
SwPhy::SetAntenna (Ptr<AntennaModel> a)
{
  NS_LOG_FUNCTION (this << a);
  m_antenna = a;
}

//----------------------------------------------------------------------
bool
SwPhy::SendPacket (Ptr<Packet> packet, bool rate)
{
  NS_LOG_FUNCTION ("");
  // RX might be interrupted by TX, but not vice versa
  if (m_state == TX) 
    {
      NS_LOG_DEBUG ("Already in transmission mode");
      return false;
    }
  
  m_state = TX;
  //Time txDuration;
  if (rate) // transmit packet with data rate
    {
      txDuration = CalTxDuration (0, packet->GetSize (), 
                        m_mac->GetBasicRate (), m_mac->GetDataRate ());
    }
  else // transmit packets (e.g. RTS, CTS) with basic rate
    {
      txDuration = CalTxDuration (packet->GetSize (), 0, 
                        m_mac->GetBasicRate (), m_mac->GetDataRate ());
    }
  
   if(m_state == TX) 
   {
      //send down
      NS_ASSERT (m_channel);
      Ptr<SwSpectrumSignalParameters> txParams = Create<SwSpectrumSignalParameters> ();
      txParams->duration = txDuration;
      txParams->txPhy = GetObject<SpectrumPhy> ();
      txParams->psd = m_txPsd;
      txParams->txAntenna = m_antenna;
      Ptr<PacketBurst> pb = CreateObject<PacketBurst> ();
      pb->AddPacket (packet);
      txParams->packetBurst = pb;
      m_channel->StartTx (txParams);
      NS_LOG_DEBUG (this << "DURATION FOR TRANSMISSION" << txParams->duration );
      Simulator::Schedule (txDuration, &SwPhy::SendPacketDone, this, packet);
      return true;
   }
  return false;
}

void 
SwPhy::SendPacketDone (Ptr<Packet> packet)
{
  NS_LOG_FUNCTION ("");
  if(m_state == COLL)
  {
        NS_LOG_UNCOND (Simulator::Now()<< " Packet transmission aborted - COLLISION (Wifi)" << " Node: " << m_mac->GetAddress());
  }
  if(txDuration < duration)
  {
        m_state = COLL;
        if(m_EndCTRXRequest.IsRunning())
        {
        m_EndCTRXRequest.Cancel();
        }
        Simulator::Schedule (duration-txDuration, &SwPhy::ChangeState, this, IDLE);
  }
  else
  {
        m_state = IDLE;
  }
  txDuration = Seconds(0);
  m_mac->SendPacketDone (packet);
}

void
SwPhy::StartRx (Ptr<SpectrumSignalParameters> spectrumRxParams)
{
  NS_LOG_FUNCTION (this << spectrumRxParams << "TIME = " << Simulator::Now());
  SwSpectrumValueHelper psdHelper;

  Ptr<SwSpectrumSignalParameters> SwRxParams = DynamicCast<SwSpectrumSignalParameters> (spectrumRxParams);
  Ptr<LrWpanSpectrumSignalParameters> LrWpanRxParams = DynamicCast<LrWpanSpectrumSignalParameters> (spectrumRxParams);

  if (SwRxParams != 0)//This is a WiFi packet.
    {       
        if(!m_EndCTRXRequest.IsRunning())//I haven't received a ZigBee packet. 
        {
                if(m_state == TX)
                {
                        NS_LOG_DEBUG ("Drop packet due to half-duplex. Already Transmitting");
                        m_state = COLL;
                }
                else 
                {
                        m_state = RX;
                }  
                if(m_EndRxRequest.IsRunning())//I have received a WiFi packet immediately before
                {      
		        NS_LOG_DEBUG (this << "*****************COLLISSION HAPPENED ********************");
                        m_state = COLL;
                        if(SwRxParams->duration > duration) 
                        {
                                m_csBusy = false;
                                m_EndRxRequest.Cancel();
                        }   
                        else
                        {
                                return; 
                        } 
                 } 
                  
                 duration = SwRxParams->duration;
              
                Ptr<Packet> p = (SwRxParams->packetBurst->GetPackets ()).front ();
                if (m_csBusy == false)
                {
                    m_csBusy = true;
                    m_pktRx = p;
                    m_mac->ReceivePacket (this, p);
                }
                m_rxPsd = SwRxParams->psd;
                m_rxTotalPower = psdHelper.TotalAvgPower (*m_rxPsd);
                m_EndRxRequest = Simulator::Schedule (duration, &SwPhy::EndRx, this);
        }
        else//I have received a ZigBee packet immediately before 
        {
	        NS_LOG_DEBUG (this << "*********COEXISTENCE COLLISSION HAPPENED ********************");
                if(SwRxParams->duration > duration)
                {
                        duration = SwRxParams->duration;
                        if(txDuration > duration)
                                duration = txDuration;
                        m_EndCTRXRequest.Cancel();
                        m_EndCTRXRequest = Simulator::Schedule (duration, &SwPhy::ChangeState, this, IDLE);
                }
                return;
        }
    }
  else//This is a ZigBee packet.
    {
      NS_LOG_DEBUG (this << "802.15.4 transmission going on, BUSY channel");
      m_state = COLL;
      
      if(m_EndRxRequest.IsRunning() || m_EndCTRXRequest.IsRunning())//I have received a WiFi/Zigbee packet immediately before
      {      
		NS_LOG_DEBUG (this << "*********COEXISTENCE COLLISSION HAPPENED ********************");
                if(LrWpanRxParams->duration > duration) 
                {
                        duration = LrWpanRxParams->duration;
                        if(txDuration > duration)
                                duration = txDuration;
                        //m_csBusy = false;
                        if(m_EndRxRequest.IsRunning())
                                m_EndRxRequest.Cancel();
                        if(m_EndCTRXRequest.IsRunning())
                                m_EndCTRXRequest.Cancel();
                }   
                else
                {
                        //Simulator::Schedule (duration, &SwPhy::ChangeState, this, IDLE);
                        return;                        
                }
       } 
       else
       {
                duration = LrWpanRxParams->duration;
                if(txDuration > duration)
                        duration = txDuration;  
       }      
       Ptr<Packet> p = (LrWpanRxParams->packetBurst->GetPackets ()).front ();
       m_mac->ReceivePacket (this, p);       
       if (m_csBusy == false)
       {
                m_csBusy = true;
                m_pktRx = p;
                m_mac->ReceivePacket (this, p);
       } 
       m_EndCTRXRequest = Simulator::Schedule (duration, &SwPhy::ChangeState, this, IDLE);
    }
  
  NS_LOG_LOGIC (this << " state: " << m_state << " m_rxTotalPower: " << m_rxTotalPower);
}

void SwPhy::EndRx ()
{
  NS_LOG_FUNCTION (this);
  m_csBusy = false;
  if (m_state != RX)
    {
      NS_LOG_INFO ("Drop packet due to state. Reception failure");
      return;
    }

  SwSpectrumValueHelper psdHelper;
 
  if (m_pktRx != 0)
    {
      double sinr = psdHelper.TotalAvgPower (*m_rxPsd) / psdHelper.TotalAvgPower (*m_noise);

      if (sinr > m_sinrTh && m_state != COLL)
        {
          NS_LOG_DEBUG ("Reception success for packet: sinr " << sinr);
          m_state = IDLE;
          m_mac->ReceivePacketDone (this, m_pktRx, true);
          m_rxTotalPower = 0;
          m_pktRx = 0;
          duration = Seconds (0);
          return;
        }
      else
        {
          /* Failure */
          NS_LOG_DEBUG ("Reception failure for packet sinr " << sinr);
        }
    }
  
  if (!m_csBusy) // set MAC state IDLE
    {
      m_state = IDLE;
      m_mac->ReceivePacketDone (this, m_pktRx, false);
    }
}

bool 
SwPhy::IsIdle ()
{
  if ((m_state == IDLE) && (m_csBusy == false)) 
  { 
        return true; 
  }
  return false;
}

void
SwPhy::SetTxPowerSpectralDensity (Ptr<SpectrumValue> txPsd)
{
  NS_LOG_FUNCTION (this << txPsd);
  NS_ASSERT (txPsd);
  m_txPsd = txPsd;
  NS_LOG_INFO ("\t computed tx_psd: " << *txPsd << "\t stored tx_psd: " << *m_txPsd);
}

void
SwPhy::SetNoisePowerSpectralDensity (Ptr<const SpectrumValue> noisePsd)
{
  NS_LOG_FUNCTION (this << noisePsd);
  NS_LOG_INFO ("\t computed noise_psd: " << *noisePsd );
  NS_ASSERT (noisePsd);
  m_noise = noisePsd;
}

Ptr<const SpectrumValue>
SwPhy::GetNoisePowerSpectralDensity (void)
{
  NS_LOG_FUNCTION (this);
  return m_noise;
}
void 
SwPhy::ChangeState (State s)
{
  NS_LOG_FUNCTION (this);
  duration = Seconds (0);
  m_state = s;
  m_csBusy = false;
  m_mac->ReceivePacketDone (this, m_pktRx, false);
}

Time
SwPhy::CalTxDuration (uint32_t basicSize, uint32_t dataSize, double basicRate, double dataRate)
{
  double_t txHdrTime = (double)(m_headerSize + basicSize + m_trailerSize) * 8.0 / basicRate;
  double_t txMpduTime = (double)dataSize * 8.0 / dataRate;
  return m_preambleDuration + Seconds (txHdrTime) + Seconds (txMpduTime);
}

} // namespace ns3
