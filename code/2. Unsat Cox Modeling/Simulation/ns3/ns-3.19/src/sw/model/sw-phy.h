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

#ifndef SW_PHY_H
#define SW_PHY_H

#include "ns3/simulator.h"
#include "ns3/event-id.h"
#include "sw-mac.h"
#include "sw-phy.h"
#include <ns3/packet-burst.h>
#include <ns3/spectrum-phy.h>
#include <ns3/spectrum-channel.h>
#include <ns3/spectrum-interference.h>
#include <ns3/spectrum-value.h>
#include <ns3/antenna-model.h>
#include "sw-spectrum-value-helper.h"
#include "sw-spectrum-signal-parameters.h"
#include "ns3/lr-wpan-spectrum-signal-parameters.h"

namespace ns3 {
/**
 * \ingroup sw
 *
 * Make SwPhy a SpectrumPhy so we can enable the eventual modeling of
 * device interference
 */

class SwPhy : public SpectrumPhy
{
public:

  enum State
    {
      IDLE, TX, RX, COLL
    };
  
  SwPhy ();
  virtual ~SwPhy ();
  void Clear ();
  void ChangeState (State s);
  static TypeId GetTypeId (void);
 
  void SetMobility (Ptr<MobilityModel> m);
  Ptr<MobilityModel> GetMobility ();
  void SetDevice (Ptr<NetDevice> device);
  Ptr<NetDevice> GetDevice ();
  void SetMac (Ptr<SwMac> mac);
  void SetChannel (Ptr<SpectrumChannel> channel);
  void SetTxPower (double dBm);

  void SetAntenna (Ptr<AntennaModel> a);
  Ptr<AntennaModel> GetRxAntenna ();

  Ptr<SpectrumChannel> GetChannel ();
  virtual Ptr<const SpectrumModel> GetRxSpectrumModel () const;
  /**
   * set the Power Spectral Density of outgoing signals in W/Hz.
   *
   * @param txPsd
   */
  void SetTxPowerSpectralDensity (Ptr<SpectrumValue> txPsd);

  /**
   * \brief set the noise power spectral density
   * @param noisePsd the Noise Power Spectral Density in power units
   * (Watt, Pascal...) per Hz.
   */
  void SetNoisePowerSpectralDensity (Ptr<const SpectrumValue> noisePsd);

  /**
   * \brief get the noise power spectral density
   * @return the Noise Power Spectral Density
   */
  Ptr<const SpectrumValue> GetNoisePowerSpectralDensity (void);

  /**
    * Notify the SpectrumPhy instance of an incoming waveform
    *
    * @param params the SpectrumSignalParameters associated with the incoming waveform
    */
  virtual void StartRx (Ptr<SpectrumSignalParameters> params);

  Mac48Address GetAddress ();
  double GetRxPowerTh ();
  double GetTxPower ();
  
  bool SendPacket (Ptr<Packet> packet, bool rate);
  void SendPacketDone (Ptr<Packet> packet);
  void EndRx ();

  bool IsIdle ();
  Time CalTxDuration (uint32_t basicSize, uint32_t dataSize, double basicRate, double dataRate);

private:
  State m_state;
  Ptr<AntennaModel> m_antenna;
  Ptr<MobilityModel> m_mobility;
  Ptr<NetDevice> m_device;
  Ptr<SwMac> m_mac;
  Ptr<SpectrumChannel> m_channel;
  Ptr<SpectrumValue> m_txPsd;
  Ptr<const SpectrumValue> m_rxPsd;
  Ptr<const SpectrumValue> m_noise;
  
  Time txDuration;
  Time duration;
  double m_rxTotalPower;
  Ptr<Packet> m_pktRx;
  Time m_preambleDuration;
  uint32_t m_trailerSize;
  uint32_t m_headerSize;
  
  double m_txPower;  // transmission power (dBm)
  double m_sinrTh;   // SINR threshold
  double m_csTh;     // carrier sense threshold (dBm)
  bool m_csBusy;
  Time m_csBusyEnd;
  EventId m_EndRxRequest;
  EventId m_EndCTRXRequest;

protected:
};

} // namespace ns3

#endif // SW_PHY_H
