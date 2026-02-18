/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2014-2026 Biswajit Banerjee, Parresia Research Limited, NZ
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

#include <sci_defs/bits_defs.h>
#include <sci_defs/malloc_defs.h>

#include <Core/Malloc/AllocOS.h>
#include <Core/Malloc/AllocPriv.h>

#ifdef __APPLE__
#include <sys/types.h>
#endif

#ifndef _WIN32
#include <sys/mman.h>
#include <unistd.h>
#endif

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>

namespace Uintah {

#ifndef DISABLE_SCI_MALLOC

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
inline constexpr int ALIGN = 32;

// On Linux and Solaris, mmap() expects a char* rather than void*.
#if defined(sun) || defined(__linux)
using MmapType = char;
#else
using MmapType = void;
#endif

// ---------------------------------------------------------------------------
// Module-level state
// ---------------------------------------------------------------------------
static int devzero_fd = -1;

// ---------------------------------------------------------------------------
// OSHunk::alloc
// ---------------------------------------------------------------------------
[[nodiscard]] OSHunk*
OSHunk::alloc(std::size_t size, bool returnable, Allocator* allocator)
{
  // Compute padding needed so that hunk->data is ALIGN-byte aligned.
  const unsigned long hunk_offset =
    (sizeof(OSHunk) % ALIGN != 0) ? (ALIGN - sizeof(OSHunk) % ALIGN) : 0UL;

  const std::size_t asize = size + sizeof(OSHunk) + hunk_offset;

  void* ptr = nullptr;

  if (returnable) {
    // Use mmap for returnable (big-object) hunks.
    if (devzero_fd == -1) {
      devzero_fd = open("/dev/zero", O_RDWR);
      if (devzero_fd == -1) {
        std::fprintf(stderr,
                     "Error opening /dev/zero: errno=%d\n", errno);
        std::abort();
      }
    }
    ptr = mmap(nullptr, asize,
               PROT_READ | PROT_WRITE, MAP_PRIVATE,
               devzero_fd, 0);
  } else {
    // Use sbrk for non-returnable (pool) hunks, ensuring alignment first.
    void* current        = sbrk(0);
    unsigned long offset = reinterpret_cast<unsigned long>(current) % ALIGN;
    if (offset != 0) {
      sbrk(static_cast<long>(ALIGN - offset));
    }
    ptr = sbrk(static_cast<long>(asize));
  }

  // Both mmap and sbrk signal failure by returning (void*)-1.
  if (ptr == reinterpret_cast<void*>(-1L)) {
    std::fprintf(stderr,
                 "Error allocating memory (%zu bytes requested)\nmmap: errno=%d\n",
                 asize, errno);

    if (allocator) {
      std::fprintf(stderr,
                   "Allocator was using %zu bytes.\n", allocator->sizealloc);

      // If the allocator is already dying, exit immediately to avoid
      // an infinite loop.
      if (allocator->dieing) {
        std::exit(1);
      }

      // Mark dying and release the lock so that any remaining
      // allocations during shutdown have a chance to succeed.
      allocator->dieing = true;
      allocator->noninline_unlock();
    }
    std::abort();
  }

  auto* hunk = static_cast<OSHunk*>(ptr);

  // Place data immediately after the OSHunk header, with alignment padding.
  hunk->data = static_cast<void*>(
    reinterpret_cast<char*>(hunk + 1) + hunk_offset);

  hunk->next       = nullptr;
  hunk->ninuse     = 0;
  hunk->len        = size;
  hunk->alloc_len  = asize;
  hunk->returnable = returnable;

  return hunk;
}

#ifdef MALLOC_TRACE
#include <MallocTraceOff.h>
#endif

// ---------------------------------------------------------------------------
// OSHunk::free
// ---------------------------------------------------------------------------
void
OSHunk::free(OSHunk* hunk)
{
  if (!hunk->returnable) {
    std::fprintf(stderr,
                 "Attempt to return a non-returnable memory hunk!\n");
    std::abort();
  }

  const std::size_t len = hunk->alloc_len;

  // munmap can spuriously fail; retry up to 10 times before giving up.
  for (int attempt = 0; attempt < 10; ++attempt) {
    if (munmap(static_cast<MmapType*>(static_cast<void*>(hunk)), len) == 0) {
      return;
    }
  }

  std::fprintf(stderr,
               "Error unmapping memory\nmunmap: errno=%d\n", errno);
  std::fprintf(stderr, "Unmap failed - leaking memory\n");
}

#ifdef MALLOC_TRACE
#include <MallocTraceOn.h>
#endif

#endif // !DISABLE_SCI_MALLOC

} // namespace Uintah