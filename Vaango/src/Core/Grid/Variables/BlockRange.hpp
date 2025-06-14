/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef UINTAH_HOMEBREW_BLOCK_RANGE_HPP
#define UINTAH_HOMEBREW_BLOCK_RANGE_HPP

#include <sci_defs/kokkos_defs.h>

#ifdef UINTAH_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // UINTAH_ENABLE_KOKKOS

#include <cstddef>
#include <vector>

namespace Uintah {

// #if defined( UINTAH_ENABLE_KOKKOS )
// typedef Kokkos::View<IntVector*> BC_List;
// #else
template<typename myIntVector>
using BC_List = std::vector<myIntVector>&;
// #endif

class BlockRange
{
public:
  enum
  {
    rank = 3
  };

  template<typename ArrayType>
  BlockRange(ArrayType const& c0, ArrayType const& c1)
  {
    for (int i = 0; i < rank; ++i) {
      m_offset[i] = c0[i] < c1[i] ? c0[i] : c1[i];
      m_dim[i]    = (c0[i] < c1[i] ? c1[i] : c0[i]) - m_offset[i];
    }
  }

  int
  begin(int r) const
  {
    return m_offset[r];
  }
  int
  end(int r) const
  {
    return m_offset[r] + m_dim[r];
  }

  size_t
  size() const
  {
    size_t result = 1u;
    for (int i = 0; i < rank; ++i) {
      result *= m_dim[i];
    }
    return result;
  }

private:
  int m_offset[rank];
  int m_dim[rank];
};

template<typename Functor>
void
serial_for(BlockRange const& r, const Functor& f)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  for (int k = kb; k < ke; ++k) {
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        f(i, j, k);
      }
    }
  }
}

template<typename myIntVector, typename Functor>
inline void
parallel_for(BC_List<myIntVector> iterSpace,
             const unsigned int bc_size,
             const Functor& functor)
// parallel_for( std::vector<IntVector> iterSpace)
{
  for (unsigned int iblock = 0; iblock < bc_size; ++iblock) {
    const int i = iterSpace[iblock][0];
    const int j = iterSpace[iblock][1];
    const int k = iterSpace[iblock][2];
    functor(i, j, k);
  }
}

#if defined(UINTAH_ENABLE_KOKKOS)

template<typename Functor>
void
parallel_for(BlockRange const& r, const Functor& f)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  Kokkos::parallel_for(
    Kokkos::RangePolicy<Kokkos::OpenMP, int>(kb, ke).set_chunk_size(2),
    KOKKOS_LAMBDA(int k) {
      for (int j = jb; j < je; ++j) {
        for (int i = ib; i < ie; ++i) {
          f(i, j, k);
        }
      }
    });
};

template<typename Functor, typename Option>
void
parallel_for(BlockRange const& r, const Functor& f, const Option& op)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  Kokkos::parallel_for(
    Kokkos::RangePolicy<Kokkos::OpenMP, int>(kb, ke).set_chunk_size(2),
    KOKKOS_LAMBDA(int k) {
      for (int j = jb; j < je; ++j) {
        for (int i = ib; i < ie; ++i) {
          f(op, i, j, k);
        }
      }
    });
};

template<typename Functor, typename ReductionType>
void
parallel_reduce_sum(BlockRange const& r, const Functor& f, ReductionType& red)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  ReductionType tmp = red;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<Kokkos::OpenMP, int>(kb, ke).set_chunk_size(2),
    KOKKOS_LAMBDA(int k, ReductionType& tmp) {
      for (int j = jb; j < je; ++j) {
        for (int i = ib; i < ie; ++i) {
          f(i, j, k, tmp);
        }
      }
    },
    tmp);
  red = tmp;
};

template<typename Functor, typename ReductionType>
void
parallel_reduce_min(BlockRange const& r, const Functor& f, ReductionType& red)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  ReductionType tmp = red;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<Kokkos::OpenMP, int>(kb, ke).set_chunk_size(2),
    KOKKOS_LAMBDA(int k, ReductionType& tmp) {
      for (int j = jb; j < je; ++j) {
        for (int i = ib; i < ie; ++i) {
          f(i, j, k, tmp);
        }
      }
    },
    Kokkos::Experimental::Min<ReductionType>(tmp));
  red = tmp;
};

#else

template<typename Functor>
void
parallel_for(BlockRange const& r, const Functor& f)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  for (int k = kb; k < ke; ++k) {
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        f(i, j, k);
      }
    }
  }
};

template<typename Functor, typename Option>
void
parallel_for(BlockRange const& r, const Functor& f, const Option& op)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  for (int k = kb; k < ke; ++k) {
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        f(op, i, j, k);
      }
    }
  }
};

template<typename Functor, typename ReductionType>
void
parallel_reduce_sum(BlockRange const& r, const Functor& f, ReductionType& red)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  ReductionType tmp = red;
  for (int k = kb; k < ke; ++k) {
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        f(i, j, k, tmp);
      }
    }
  }
  red = tmp;
};

template<typename Functor, typename ReductionType>
void
parallel_reduce_min(BlockRange const& r, const Functor& f, ReductionType& red)
{
  const int ib = r.begin(0);
  const int ie = r.end(0);
  const int jb = r.begin(1);
  const int je = r.end(1);
  const int kb = r.begin(2);
  const int ke = r.end(2);

  ReductionType tmp = red;
  for (int k = kb; k < ke; ++k) {
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        f(i, j, k, tmp);
      }
    }
  }
  red = tmp;
};

#endif

} // namespace Uintah

#endif // UINTAH_HOMEBREW_BLOCK_RANGE_HPP
