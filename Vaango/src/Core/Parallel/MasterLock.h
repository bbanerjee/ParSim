/*
 * The MIT License
 *
 * Copyright (c) 1997-2022 The University of Utah
 * Copyright (c) 2018-2023 Parresia Research Limited, NZ
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

#ifndef __CORE_PARALLEL_MASTERLOCK_H__
#define __CORE_PARALLEL_MASTERLOCK_H__

#include <sci_defs/kokkos_defs.h>

#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Uintah {

class MasterLock
{

  // Specific to OpenMP- and std::thread-based implementations

  // This lock can be used as follows:
  //
  // 1.) Same functionality as std::mutex, e.g., create, lock(), and unlock()
  //
  //      OR
  //
  // 2.) std::unique_lock<MasterLock>
  //
  //       std::unique_lock is a general-purpose mutex ownership wrapper
  //       allowing deferred locking, time-constrained attempts at locking,
  //       recursive locking, transfer of lock ownership, and use with condition
  //       variables.
  //
  //      OR
  //
  // 3.) std::lock_guard<MasterLock>
  //
  //       std::lock_guard is a mutex wrapper that provides a convenient
  //       RAII-style mechanism for owning a mutex for the duration of a scoped
  //       block.

public:
#if defined(_OPENMP) && defined(UINTAH_ENABLE_KOKKOS)

  // per OMP standard, a flush region without a list is implied for
  // omp_{set/unset}_lock
  void lock() { omp_set_lock(&m_lock); }
  void unlock() { omp_unset_lock(&m_lock); }

  MasterLock() { omp_init_lock(&m_lock); }
  ~MasterLock() { omp_destroy_lock(&m_lock); }

#else

  void lock() { m_mutex.lock(); }
  void unlock() { m_mutex.unlock(); }

  MasterLock() {}
  ~MasterLock() {}

#endif // UINTAH_ENABLE_KOKKOS

private:
  // disable copy, assignment, and move
  MasterLock(const MasterLock&) = delete;
  MasterLock& operator=(const MasterLock&) = delete;
  MasterLock(MasterLock&&) = delete;
  MasterLock& operator=(MasterLock&&) = delete;

#if defined(_OPENMP) && defined(UINTAH_ENABLE_KOKKOS)
  omp_lock_t m_lock;
#else
  std::mutex m_mutex;
#endif // UINTAH_ENABLE_KOKKOS
};
} // end namespace Uintah

#endif // end __CORE_PARALLEL_MASTERLOCK_H__
