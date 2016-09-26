/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2006 INRIA
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
 * Author: Mathieu Lacage <mathieu.lacage@sophia.inria.fr>
 */
#include "full-ssid.h"
#include "ns3/assert.h"

namespace ns3 {

FullSsid::FullSsid ()
{
  m_length = 0;
  for (uint8_t i = 0; i < 33; i++)
    {
      m_ssid[i] = 0;
    }
}
FullSsid::FullSsid (std::string s)
{
  NS_ASSERT (s.size () < 32);
  const char *ssid = s.c_str ();
  uint8_t len = 0;
  while (*ssid != 0 && len < 32)
    {
      m_ssid[len] = *ssid;
      ssid++;
      len++;
    }
  NS_ASSERT (len <= 32);
  m_length = len;
  while (len < 33)
    {
      m_ssid[len] = 0;
      len++;
    }
}
FullSsid::FullSsid (char const ssid[32], uint8_t length)
{
  NS_ASSERT (length <= 32);
  uint8_t len = 0;
  while (len < length)
    {
      m_ssid[len] = ssid[len];
      len++;
    }
  m_length = length;
  while (len < 33)
    {
      m_ssid[len] = 0;
      len++;
    }
}
bool
FullSsid::IsEqual (const FullSsid& o) const
{
  uint8_t i = 0;
  while (i < 32
         && m_ssid[i] == o.m_ssid[i]
         && m_ssid[i] != 0)
    {
      i++;
    }
  if (m_ssid[i] != o.m_ssid[i])
    {
      return false;
    }
  return true;
}
bool
FullSsid::IsBroadcast (void) const
{
  if (m_ssid[0] == 0)
    {
      return true;
    }
  return false;
}

char *
FullSsid::PeekString (void) const
{
  // it is safe to return a pointer to the buffer because it is
  // guaranteed to be zero-terminated.
  return (char *)m_ssid;
}

FullWifiInformationElementId
FullSsid::ElementId () const
{
  return IE_SSID;
}

uint8_t
FullSsid::GetInformationFieldSize () const
{
  return m_length;
}

void
FullSsid::SerializeInformationField (Buffer::Iterator start) const
{
  NS_ASSERT (m_length <= 32);
  start.Write (m_ssid, m_length);
}
uint8_t
FullSsid::DeserializeInformationField (Buffer::Iterator start,
                                   uint8_t length)
{
  m_length = length;
  NS_ASSERT (m_length <= 32);
  start.Read (m_ssid, m_length);
  return length;
}

ATTRIBUTE_HELPER_CPP (FullSsid);

std::ostream &
operator << (std::ostream &os, const FullSsid &ssid)
{
  os << ssid.PeekString ();
  return os;
}

std::istream &operator >> (std::istream &is, FullSsid &ssid)
{
  std::string str;
  is >> str;
  ssid = FullSsid (str.c_str ());
  return is;
}


} // namespace ns3
