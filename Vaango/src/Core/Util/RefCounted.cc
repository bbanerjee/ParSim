/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <Core/Util/RefCounted.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/MasterLock.h>

#include <Core/Util/Assert.h>
#include <Core/Util/FancyAssert.h>

#include <atomic>

namespace Uintah {

static const int NLOCKS = 1024;
static bool initialized = false;

static Uintah::MasterLock* locks[NLOCKS];
static Uintah::MasterLock initlock{};

static std::atomic<long long> nextIndex{0};
static std::atomic<long long> freeIndex{0};

RefCounted::RefCounted()
  : d_refCount(0)
{
  if (!initialized) {
    initlock.lock();
    {
      for (int i = 0; i < NLOCKS; i++) {
        locks[i] = scinew Uintah::MasterLock();
      }
      initialized = true;
    }
    initlock.unlock();
  }
  d_lockIndex = (nextIndex.fetch_add(1, std::memory_order_seq_cst)) % NLOCKS;
  ASSERT(d_lockIndex >= 0);
}

/*
RefCounted::~RefCounted() noexcept(false)
{
  ASSERTEQ(d_refCount, 0);
  int index = freeIndex.load(std::memory_order_seq_cst);

  if (index == nextIndex.load(std::memory_order_seq_cst)) {
    initlock.lock();
    {
      if (freeIndex.store(nextIndex.load(std::memory_order_seq_cst)),
                          std::memory_order_seq_cst) {
        initialized = false;
        for (int i = 0; i < NLOCKS; i++) {
          delete locks[i];
          locks[i] = nullptr;
        }

      }
    }
    initlock.unlock();
  }
}
*/

// nextIndex : tracks the total number of RefCounted objects created
// freeIndex : appears to track objects that have been destroyed
// The destructor should increment freeIndex and check if all objects are destroyed
// The issue is that you're trying to use store() in an if condition, 
// but you should be using fetch_add() to increment freeIndex and 
// then check if it equals nextIndex.
/*
RefCounted::~RefCounted() noexcept(false)
{
  ASSERTEQ(d_refCount, 0);

  // Increment freeIndex and get the new value
  long long newFreeIndex = freeIndex.fetch_add(1, std::memory_order_seq_cst) + 1;
  
  // Check if all objects have been destroyed
  if (newFreeIndex == nextIndex.load(std::memory_order_seq_cst)) {
    initlock.lock();
    {
      // Double-check inside the lock to avoid race conditions
      if (freeIndex.load(std::memory_order_seq_cst) == nextIndex.load(std::memory_order_seq_cst)) {
        initialized = false;
        for (int i = 0; i < NLOCKS; i++) {
          delete locks[i];
          locks[i] = nullptr;
        }
        // Reset counters for potential reuse
        nextIndex.store(0, std::memory_order_seq_cst);
        freeIndex.store(0, std::memory_order_seq_cst);
      }
    }
    initlock.unlock();
  }
}
*/

RefCounted::~RefCounted() noexcept(false)
{
  ASSERTEQ(d_refCount, 0);
  
  // Increment freeIndex
  long long newFreeIndex = freeIndex.fetch_add(1, std::memory_order_seq_cst) + 1;
  long long currentNext = nextIndex.load(std::memory_order_seq_cst);
  
  // Check if this is the last object being destroyed
  if (newFreeIndex == currentNext) {
    initlock.lock();
    {
      // Use compare_exchange to ensure atomicity
      long long expected = currentNext;
      if (freeIndex.compare_exchange_strong(expected, currentNext, std::memory_order_seq_cst) &&
          expected == currentNext) {
        initialized = false;
        for (int i = 0; i < NLOCKS; i++) {
          delete locks[i];
          locks[i] = nullptr;
        }
        // Reset for potential reuse
        nextIndex.store(0, std::memory_order_seq_cst);
        freeIndex.store(0, std::memory_order_seq_cst);
      }
    }
    initlock.unlock();
  }
}

/*
// If lock already exists
RefCounted::~RefCounted() noexcept(false)
{
  ASSERTEQ(d_refCount, 0);
  int index = freeIndex.load(std::memory_order_seq_cst);
  
  if (index == nextIndex.load(std::memory_order_seq_cst)) {
    initlock.lock();
    {
      // Since we're locked, we can safely check and update
      if (freeIndex.load(std::memory_order_seq_cst) == nextIndex.load(std::memory_order_seq_cst)) {
        initialized = false;
        for (int i = 0; i < NLOCKS; i++) {
          delete locks[i];
          locks[i] = nullptr;
        }
        // Only reset freeIndex after cleanup is done
        freeIndex.store(0, std::memory_order_seq_cst);  // or whatever reset value you need
      }
    }
    initlock.unlock();
  }
}
*/

void
RefCounted::addReference() const
{
  locks[d_lockIndex]->lock();
  {
    d_refCount++;
  }
  locks[d_lockIndex]->unlock();
}

bool
RefCounted::removeReference() const
{
  bool status;
  locks[d_lockIndex]->lock();
  {
    status = (--d_refCount == 0);
    ASSERT(d_refCount >= 0);
  }
  locks[d_lockIndex]->unlock();
  return status;
}

} // end namespace Uintah