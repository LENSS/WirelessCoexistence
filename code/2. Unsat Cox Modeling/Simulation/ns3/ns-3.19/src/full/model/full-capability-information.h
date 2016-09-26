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
#ifndef FULL_CAPABILITY_INFORMATION_H
#define FULL_CAPABILITY_INFORMATION_H

#include <stdint.h>
#include "ns3/buffer.h"

namespace ns3 {

/**
 * \ingroup wifi
 *
 *
 */
class FullCapabilityInformation
{
public:
  FullCapabilityInformation ();

  void SetEss (void);
  void SetIbss (void);

  bool IsEss (void) const;
  bool IsIbss (void) const;

  uint32_t GetSerializedSize (void) const;
  Buffer::Iterator Serialize (Buffer::Iterator start) const;
  Buffer::Iterator Deserialize (Buffer::Iterator start);
private:
  bool Is (uint8_t n) const;
  void Set (uint8_t n);
  void Clear (uint8_t n);
  uint16_t m_capability;
};

} // namespace ns3

#endif /* FULL_CAPABILITY_INFORMATION_H */