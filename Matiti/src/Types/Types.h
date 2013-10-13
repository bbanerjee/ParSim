#ifndef MATITI_TYPES_H
#define MATITI_TYPES_H

#include <cstdint>
#include <array>

namespace Matiti {
  
  typedef uint8_t u8;
  typedef uint32_t u32;
  typedef uint64_t u64;
  typedef uint64_t u128;
  typedef int64_t long64;

  typedef std::array<int, 3> IntArray3;
  typedef std::array<double, 3> Array3;
}

#endif
