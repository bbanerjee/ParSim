#ifndef PERI_ELEMENT_H
#define PERI_ELEMENT_H

#include <vector>

namespace pd {

  struct PeriElement
  {
    std::vector<int> nodes = std::vector<int>(8, 0);
    const int& operator[](int idx) const
    {
      return (idx < 0 || idx > 7) ? nodes[0] : nodes[idx];
    }
    int& operator[](int idx)
    {
      return (idx < 0 || idx > 7) ? nodes[0] : nodes[idx];
    }
  };

} // end pd

#endif
