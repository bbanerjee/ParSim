#ifndef PERI_ELEMENT_H
#define PERI_ELEMENT_H

#include <Core/Types/IntegerTypes.h>
#include <vector>

namespace pd {

  struct PeriElement
  {
    std::vector<ParticleID> nodes = std::vector<ParticleID>(8, 0);
    const ParticleID& operator[](int idx) const
    {
      return (idx < 0 || idx > 7) ? nodes[0] : nodes[idx];
    }
    ParticleID& operator[](int idx)
    {
      return (idx < 0 || idx > 7) ? nodes[0] : nodes[idx];
    }
  };

} // end pd

#endif
