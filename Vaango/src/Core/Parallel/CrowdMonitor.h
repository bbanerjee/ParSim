/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2021-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef __CORE_PARALLEL_CROWDMONITOR_H__
#define __CORE_PARALLEL_CROWDMONITOR_H__

#include <Core/Parallel/MasterLock.h>

#include <atomic>
#include <thread>

namespace Uintah {

template <typename Tag = void>
class CrowdMonitor
{

public:

  enum MonitorType { READER, WRITER };

  CrowdMonitor( MonitorType t)
    : m_type(t)
  {
    if (m_type == READER) {
      std::unique_lock<Uintah::MasterLock>{ s_mutex };
      s_count.fetch_add(1, std::memory_order_relaxed);
    }
    else {
      s_mutex.lock();
//       spin_wait( s_count.fetch( std::memory_order_relaxed) != 0 );
    }
  }

  ~CrowdMonitor()
  {
    if ( m_type == READER) {
      s_count.fetch_add(-1, std::memory_order_relaxed);
    }
    else {
      s_mutex.unlock();
    }
  }

  // disable copy, assignment, and move
  CrowdMonitor( const CrowdMonitor & )            = delete;
  CrowdMonitor& operator=( const CrowdMonitor & ) = delete;
  CrowdMonitor( CrowdMonitor && )                 = delete;
  CrowdMonitor& operator=( CrowdMonitor && )      = delete;

private:

  static Uintah::MasterLock s_mutex;
  static std::atomic<int>   s_count;
  MonitorType               m_type;

};

template <typename Tag> Uintah::MasterLock CrowdMonitor<Tag>::s_mutex = {};
template <typename Tag> std::atomic<int>   CrowdMonitor<Tag>::s_count = {0};

} // end namespace Uintah

#endif // end __CORE_PARALLEL_CROWDMONITOR_H__
