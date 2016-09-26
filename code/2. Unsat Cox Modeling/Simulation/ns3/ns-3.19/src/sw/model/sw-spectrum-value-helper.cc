/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011 The Boeing Company
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Gary Pei <guangyu.pei@boeing.com>
 */
#include "ns3/log.h"

#include "sw-spectrum-value-helper.h"
NS_LOG_COMPONENT_DEFINE ("SwSpectrumValueHelper");

namespace ns3 {

Ptr<SpectrumModel> g_SwSpectrumModel;

class SwSpectrumModelInitializer
{
public:
  SwSpectrumModelInitializer ()
  {
    NS_LOG_FUNCTION (this);

    Bands bands;
    // 1 MHz resolution from 2400 to 2483 MHz (center freq)
    for (int i = -1; i < 85; i++)
      {
        BandInfo bi;
        bi.fl = 2400.5e6 + i * 1.0e6;
        bi.fh = 2400.5e6 + (i + 1) * 1.0e6;
        bi.fc = (bi.fl +  bi.fh) / 2;
        bands.push_back (bi);
      }
    g_SwSpectrumModel = Create<SpectrumModel> (bands);
  }

} g_SwSpectrumModelInitializerInstance;

SwSpectrumValueHelper::SwSpectrumValueHelper ()
{
  NS_LOG_FUNCTION (this);
  m_noiseFactor = 1.0;
}

SwSpectrumValueHelper::~SwSpectrumValueHelper ()
{
}

Ptr<SpectrumValue>
SwSpectrumValueHelper::CreateTxPowerSpectralDensity (double txPower, uint32_t channel)
{
  NS_LOG_FUNCTION (this);
  Ptr<SpectrumValue> txPsd = Create <SpectrumValue> (g_SwSpectrumModel);

  // txPower is expressed in dBm. We must convert it into natural unit.
  txPower = pow (10., (txPower - 30) / 10);

  double txPowerDensity = txPower / 2.0e6;

  NS_ASSERT_MSG ((channel >= 1 && channel <= 14), "Invalid channel numbers");

  // The channel assignment is in section 6.1.2.1
  // Table 25 specify the Tx PSD
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 - 5] = txPowerDensity * 0.01; // -20 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 - 4] = txPowerDensity * 0.01; // -20 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 - 3] = txPowerDensity * 0.1; // -10 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 - 2] = txPowerDensity * 0.1; // -10 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 - 1] = txPowerDensity * 0.5; // -3 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400] = txPowerDensity; // center
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 + 1 ] = txPowerDensity * 0.5; // -3 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 + 2 ] = txPowerDensity * 0.1; // -10 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 + 3 ] = txPowerDensity * 0.1; // -10 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 + 4 ] = txPowerDensity * 0.01; // -20 dB
  (*txPsd)[2405 + 5 * (channel - 1) - 2400 + 5 ] = txPowerDensity * 0.01; // -20 dB

  return txPsd;
}

Ptr<SpectrumValue>
SwSpectrumValueHelper::CreateNoisePowerSpectralDensity (uint32_t channel)
{
  NS_LOG_FUNCTION (this);
  Ptr<SpectrumValue> noisePsd = Create <SpectrumValue> (g_SwSpectrumModel);

  static const double BOLTZMANN = 1.3803e-23;
  // Nt  is the power of thermal noise in W
  double Nt = BOLTZMANN * 290.0;
  // noise Floor (W) which accounts for thermal noise and non-idealities of the receiver
  double noisePowerDensity = m_noiseFactor * Nt;

  NS_ASSERT_MSG ((channel >= 1 && channel <= 14), "Invalid channel numbers");

  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 - 1] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 - 2] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 - 3] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 - 4] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 - 5] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 + 1] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 + 2] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 + 3] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 + 4] = noisePowerDensity;
  (*noisePsd)[2405 + 5 * (channel - 1) - 2400 + 5] = noisePowerDensity;

  return noisePsd;
}

double
SwSpectrumValueHelper::TotalAvgPower (const SpectrumValue &psd)
{
  NS_LOG_FUNCTION (this);
  double totalAvgPower = 0.0;

  // numerically integrate to get area under psd using
  // 1 MHz resolution from 2400 to 2483 MHz (center freq)

  totalAvgPower = Sum (psd * 1.0e6);
  return totalAvgPower;
}

} // namespace ns3
