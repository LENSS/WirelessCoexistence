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
#ifndef SW_SPECTRUM_VALUE_HELPER_H
#define SW_SPECTRUM_VALUE_HELPER_H

#include <ns3/spectrum-value.h>
#include <cmath>

namespace ns3 {

/**
 * \ingroup Sw
 *
 * \brief This class defines all functions to create spectrum model for Sw
 */
class SwSpectrumValueHelper
{
public:
  SwSpectrumValueHelper ();
  virtual ~SwSpectrumValueHelper ();

  /**
   * \brief create spectrum value
   * \param txPower the power transmission in dBm
   * \param channel the channel number per IEEE802.11
   * \return a Ptr to a newly created SpectrumValue instance
   */
  Ptr<SpectrumValue> CreateTxPowerSpectralDensity (double txPower, uint32_t channel);

  /**
   * \brief create spectrum value for noise
   * \param channel the channel number per IEEE802.11
   * \return a Ptr to a newly created SpectrumValue instance
   */
  Ptr<SpectrumValue> CreateNoisePowerSpectralDensity (uint32_t channel);

  /**
   * \brief total average power of the signal is the integral of the PSD
   * \param power spectral density
   * \return total power (using composite trap. rule to numerally integrate
   */
  double TotalAvgPower (const SpectrumValue &psd);

private:
  double m_noiseFactor;

};

} //namespace ns3

#endif /*  SW_SPECTRUM_VALUE_HELPER_H */
