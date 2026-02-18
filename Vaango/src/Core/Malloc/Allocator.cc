/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2014-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifdef __INTEL_COMPILER
// Disable the fprintf warning that appears everywhere in this file on icpc.
#pragma warning(disable : 181)
#endif

/* TODO:
6) Destroy allocators
*/

#include <sci_defs/bits_defs.h>
#include <Core/Malloc/Allocator.h>

// USE_LENNY_HACK: See Allocator.h for more information.

#if !defined(DISABLE_SCI_MALLOC)

#include <Core/Malloc/AllocOS.h>
#include <Core/Malloc/AllocPriv.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#ifdef SCI_PTHREAD
#include <pthread.h>
#endif

#ifndef _WIN32
#include <sys/param.h>
#include <strings.h>
#endif

namespace Uintah {

// ---------------------------------------------------------------------------
// Alignment
// ---------------------------------------------------------------------------
inline constexpr int ALIGN = 32;

// ---------------------------------------------------------------------------
// Object state constants
// ---------------------------------------------------------------------------
inline constexpr int OBJFREE            = 1;
inline constexpr int OBJINUSE           = 2;
inline constexpr int OBJFREEING         = 3;
inline constexpr int OBJMEMALIGNFREEING = 4;

// ---------------------------------------------------------------------------
// Sentinel values
// ---------------------------------------------------------------------------
inline constexpr unsigned int SENT_VAL_FREE  = 0xdeadbeef;
inline constexpr unsigned int SENT_VAL_INUSE = 0xbeefface;

// Fill pattern used to detect use-after-free / out-of-bounds writes
inline constexpr unsigned int FILL_PATTERN = 0xffff5a5au;

// ---------------------------------------------------------------------------
// Bin sizing: small objects  (granularity: 16 bytes, max: ~504 bytes)
// ---------------------------------------------------------------------------
inline constexpr std::size_t SMALL_THRESHOLD = 512 - 8;

[[nodiscard]] constexpr std::size_t small_bin(std::size_t size) noexcept
{
  return (size + 7) >> 4;
}
[[nodiscard]] constexpr std::size_t small_binsize(std::size_t bin) noexcept
{
  return (bin << 4) + 8;
}

inline constexpr std::size_t NSMALL_BINS = (SMALL_THRESHOLD + 8) >> 4;

// ---------------------------------------------------------------------------
// Bin sizing: medium objects  (granularity: 2 KiB, max: ~512 KiB)
// ---------------------------------------------------------------------------
inline constexpr std::size_t MEDIUM_THRESHOLD = 65536 * 8;

[[nodiscard]] constexpr std::size_t medium_bin(std::size_t size) noexcept
{
  return (size - 1) >> 11;
}
[[nodiscard]] constexpr std::size_t medium_binsize(std::size_t bin) noexcept
{
  return (bin << 11) + 2048;
}

inline constexpr std::size_t NMEDIUM_BINS = MEDIUM_THRESHOLD >> 11;

// ---------------------------------------------------------------------------
// Other allocator constants
// ---------------------------------------------------------------------------
inline constexpr std::size_t OVERHEAD            = sizeof(Tag) + 2 * sizeof(Sentinel);
inline constexpr std::size_t SMALLEST_ALLOCSIZE  = 8 * 1024;
inline constexpr std::size_t NORMAL_OS_ALLOC_SIZE = 512 * 1024;
inline constexpr std::size_t MAX_ALLOCSIZE        = 1024 * 1024 * 1024;

// AIX defines STATSIZE in sys/param.h — undefine it first.
#ifdef STATSIZE
#undef STATSIZE
#endif
inline constexpr int STATSIZE = 4096 + BUFSIZ;

// ---------------------------------------------------------------------------
// Module-level state
// ---------------------------------------------------------------------------
static char trace_buffer[STATSIZE];

Allocator* default_allocator = nullptr;

static bool do_shutdown         = false;
static int  mallocStatsAppendNum = -1;

// ---------------------------------------------------------------------------
// Helper: initialise memory when INITIALIZE_MEMORY is defined
// ---------------------------------------------------------------------------
static inline void initialize_memory(void* mem, std::size_t size) noexcept
{
#ifdef INITIALIZE_MEMORY
  std::fill_n(static_cast<unsigned char*>(mem), size, MEMORY_INIT_NUMBER);
#else
  (void)mem;
  (void)size;
#endif
}

// ---------------------------------------------------------------------------

void
AllocatorMallocStatsAppendNumber(int num)
{
  mallocStatsAppendNum = num;
}

inline std::size_t
Allocator::obj_maxsize(Tag* t) noexcept
{
  if (!t->bin) {
    return t->reqsize;
  }
  return (t->bin == &big_bin) ? (t->hunk->len - OVERHEAD) : t->bin->maxsize;
}

// ---------------------------------------------------------------------------
// Internal accounting helper
// ---------------------------------------------------------------------------
static void
account_bin(Allocator* a,
            AllocBin*  bin,
            FILE*      out,
            std::size_t& bytes_overhead,
            std::size_t& bytes_free,
            std::size_t& bytes_fragmented,
            std::size_t& bytes_inuse)
{
  for (Tag* p = bin->free; p != nullptr; p = p->next) {
    bytes_overhead += OVERHEAD;
    bytes_free     += a->obj_maxsize(p);
  }
  for (Tag* p = bin->inuse; p != nullptr; p = p->next) {
    bytes_overhead    += OVERHEAD;
    bytes_inuse       += p->reqsize;
    bytes_fragmented  += a->obj_maxsize(p) - p->reqsize;
    if (out) {
#ifdef USE_TAG_LINENUM
      fprintf(out,
              "%p: %zu bytes (%s:%d)\n",
              static_cast<void*>(static_cast<char*>(static_cast<void*>(p))
                                 + sizeof(Tag) + sizeof(Sentinel)),
              p->reqsize,
              p->tag,
              p->linenum);
#else
      fprintf(out,
              "%p: %zu bytes (%s)\n",
              static_cast<void*>(static_cast<char*>(static_cast<void*>(p))
                                 + sizeof(Tag) + sizeof(Sentinel)),
              p->reqsize,
              p->tag);
#endif
    }
  }
}

// ---------------------------------------------------------------------------
// Shutdown / stats dump
// ---------------------------------------------------------------------------
#if !defined(USE_LENNY_HACK)
static
#endif
void
shutdown()
{
  static char stat_buffer[STATSIZE];

  Allocator* a = DefaultAllocator();
  if (a->statsfile && !a->stats_out) {
    char filename[256];
    std::strcpy(filename, a->statsfile);
    if (mallocStatsAppendNum >= 0) {
      std::strcat(filename, ".");
      std::sprintf(filename + std::strlen(filename), "%i", mallocStatsAppendNum);
    }
    a->stats_out = std::fopen(filename, "w");
    std::setvbuf(a->stats_out, stat_buffer, _IOFBF, STATSIZE);
    if (!a->stats_out) {
      std::perror("fopen");
      std::fprintf(stderr,
                   "cannot open stats file: %s, will not print stats\n",
                   filename);
      a->stats_out = nullptr;
    }
  }

  if (!a->stats_out) {
    return;
  }

  if (do_shutdown) {
    // Called again — rewind the file.
    std::rewind(a->stats_out);
  }

  a->lock();

  std::fprintf(a->stats_out, "Unfreed objects:\n");

  std::size_t bytes_overhead = 0, bytes_free = 0,
              bytes_fragmented = 0, bytes_inuse = 0, bytes_inhunks = 0;

  for (int i = 0; i < static_cast<int>(NSMALL_BINS); i++) {
    account_bin(a, &a->small_bins[i], a->stats_out,
                bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);
  }
  for (int i = 0; i < static_cast<int>(NMEDIUM_BINS); i++) {
    account_bin(a, &a->medium_bins[i], a->stats_out,
                bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);
  }
  account_bin(a, &a->big_bin, a->stats_out,
              bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);

  for (OSHunk* hunk = a->hunks; hunk != nullptr; hunk = hunk->next) {
    bytes_overhead += sizeof(OSHunk);
    bytes_inhunks  += hunk->spaceleft;
  }
  for (Tag* p = a->big_bin.free;  p != nullptr; p = p->next) { bytes_overhead += sizeof(OSHunk); }
  for (Tag* p = a->big_bin.inuse; p != nullptr; p = p->next) { bytes_overhead += sizeof(OSHunk); }
  bytes_overhead += a->mysize;

  if (bytes_inuse == 0) {
    std::fprintf(a->stats_out, "None\n");
  }

  std::fprintf(a->stats_out, "statistics:\n");
  std::fprintf(a->stats_out, "alloc:\t\t\t%zu calls\n",  a->nalloc);
  std::fprintf(a->stats_out, "alloc:\t\t\t%zu bytes\n",  a->sizealloc);
  std::fprintf(a->stats_out, "free:\t\t\t%zu calls\n",   a->nfree);
  std::fprintf(a->stats_out, "free:\t\t\t%zu bytes\n",   a->sizefree);
  std::fprintf(a->stats_out, "fillbin:\t\t%zu calls\n",  a->nfillbin);
  std::fprintf(a->stats_out, "mmap:\t\t\t%zu calls\n",   a->nmmap);
  std::fprintf(a->stats_out, "mmap:\t\t\t%zu bytes\n",   a->sizemmap);
  std::fprintf(a->stats_out, "munmap:\t\t\t%zu calls\n", a->nmunmap);
  std::fprintf(a->stats_out, "munmap:\t\t\t%zu bytes\n", a->sizemunmap);
  std::fprintf(a->stats_out, "highwater alloc:\t%zu bytes\n", a->highwater_alloc);
  std::fprintf(a->stats_out, "highwater mmap:\t\t%zu bytes\n", a->highwater_mmap);
  std::fprintf(a->stats_out, "\n");
  std::fprintf(a->stats_out, "breakdown of total bytes:\n");
  std::fprintf(a->stats_out, "in use:\t\t\t%zu bytes\n",        bytes_inuse);
  std::fprintf(a->stats_out, "free:\t\t\t%zu bytes\n",          bytes_free);
  std::fprintf(a->stats_out, "fragmentation:\t\t%zu bytes\n",   bytes_fragmented);
  std::fprintf(a->stats_out, "left in mmap hunks:\t%zu\n",      bytes_inhunks);
  std::fprintf(a->stats_out, "per object overhead:\t%zu bytes\n", bytes_overhead);
  std::fprintf(a->stats_out, "\n");
  std::fprintf(a->stats_out,
               "%zu bytes missing (%zu memory objects)\n",
               a->sizealloc - a->sizefree,
               a->nalloc    - a->nfree);

  a->unlock();
  do_shutdown = true;
}

// ---------------------------------------------------------------------------
// Bin lookup
// ---------------------------------------------------------------------------
inline AllocBin*
Allocator::get_bin(std::size_t size) noexcept
{
  if (size <= SMALL_THRESHOLD) {
    return &small_bins[small_bin(size)];
  } else if (size <= MEDIUM_THRESHOLD) {
    return &medium_bins[medium_bin(size)];
  } else {
    return &big_bin;
  }
}

// ---------------------------------------------------------------------------
// Locking
// ---------------------------------------------------------------------------
#if defined(SCI_NOTHREAD) || defined(DISABLE_SCI_MALLOC)

void  Allocator::initlock() noexcept {}
void  Allocator::lock()    noexcept {}
void  Allocator::unlock()  noexcept {}
void  LockAllocator(Allocator*)   {}
void  UnLockAllocator(Allocator*) {}

#else
#ifdef SCI_PTHREAD

void
Allocator::initlock() noexcept
{
  use_rlock         = false;
  lock_count        = 0;
  owner             = 0;
  owner_initialized = false;

  static pthread_mutex_t init = PTHREAD_MUTEX_INITIALIZER;
  the_lock = init;
}

inline void
Allocator::lock() noexcept
{
  if (!use_rlock) {
    if (pthread_mutex_lock(&the_lock) != 0) {
      std::perror("Allocator::lock: pthread_mutex_lock");
      std::exit(-1);
    }
  } else {
    pthread_t me = pthread_self();
    if (owner_initialized && pthread_equal(owner, me)) {
      ++lock_count;
      return;
    }
    if (pthread_mutex_lock(&the_lock) != 0) {
      std::perror("Allocator::lock: pthread_mutex_lock");
      std::exit(-1);
    }
  }
}

inline void
Allocator::rlock() noexcept
{
  pthread_t me = pthread_self();
  if (owner_initialized && pthread_equal(owner, me)) {
    ++lock_count;
    return;
  }
  if (pthread_mutex_lock(&the_lock) != 0) {
    std::perror("Allocator::rlock: pthread_mutex_lock");
    std::exit(-1);
  }
  owner             = me;
  owner_initialized = true;
  use_rlock         = true;
  lock_count        = 1;
}

inline void
Allocator::unlock() noexcept
{
  if (!use_rlock) {
    if (pthread_mutex_unlock(&the_lock) != 0) {
      std::perror("Allocator::unlock: pthread_mutex_unlock");
      std::exit(-1);
    }
  } else {
    if (--lock_count == 0) {
      owner             = 0;
      owner_initialized = false;
      use_rlock         = false;
      if (pthread_mutex_unlock(&the_lock) != 0) {
        std::perror("Allocator::unlock: pthread_mutex_unlock");
        std::exit(-1);
      }
    }
  }
}

void LockAllocator(Allocator* a)   { a->rlock();  }
void UnLockAllocator(Allocator* a) { a->unlock(); }

#endif // SCI_PTHREAD
#endif // SCI_NOTHREAD || DISABLE_SCI_MALLOC

// ---------------------------------------------------------------------------
// Public factory / error helpers
// ---------------------------------------------------------------------------
void
MakeDefaultAllocator()
{
  if (!default_allocator) {
    default_allocator = MakeAllocator();
  }
}

void
AllocError(const char* msg)
{
  std::fprintf(stderr, "Allocator error: %s\n", msg);
  std::abort();
}

// ---------------------------------------------------------------------------
// MakeAllocator
// ---------------------------------------------------------------------------
[[nodiscard]] Allocator*
MakeAllocator()
{
#ifndef DISABLE_SCI_MALLOC
  std::size_t size = sizeof(Allocator)
                   + NSMALL_BINS  * sizeof(AllocBin)
                   + NMEDIUM_BINS * sizeof(AllocBin);

  OSHunk* alloc_hunk    = OSHunk::alloc(size, false, nullptr);
  auto*   a             = static_cast<Allocator*>(alloc_hunk->data);
  alloc_hunk->spaceleft = 0;
  alloc_hunk->next      = nullptr;
  alloc_hunk->ninuse    = 1;

  a->hunks      = alloc_hunk;
  a->small_bins = reinterpret_cast<AllocBin*>(a + 1);
  a->mysize     = size;

  for (int j = 0; j < static_cast<int>(NSMALL_BINS); j++) {
    std::size_t minsize = (j == 0) ? 0 : small_binsize(j - 1) + 1;
    a->init_bin(&a->small_bins[j], small_binsize(j), minsize);
  }

  a->medium_bins = a->small_bins + NSMALL_BINS;
  for (int i = 0; i < static_cast<int>(NMEDIUM_BINS); i++) {
    std::size_t minsize = (i == 0) ? SMALL_THRESHOLD + 1 : medium_binsize(i - 1) + 1;
    a->init_bin(&a->medium_bins[i], medium_binsize(i), minsize);
  }
  a->init_bin(&a->big_bin, MAX_ALLOCSIZE, MEDIUM_THRESHOLD + 1);

  a->strict = (getenv("MALLOC_STRICT") != nullptr) ? 1 : 0;
  a->lazy   = (getenv("MALLOC_LAZY")   != nullptr) ? 1 : 0;

  // Initialise stats
  a->nmmap          = 1;
  a->sizemmap       = size + sizeof(OSHunk);
  a->highwater_mmap = size;
  a->nalloc = a->nfree = a->sizealloc = a->sizefree = 0;
  a->nfillbin                                       = 0;
  a->nmunmap = a->sizemunmap                        = 0;
  a->highwater_alloc                                = 0;

  a->pagesize = getpagesize();
  a->initlock();

  bool atexit_added = false;

  // --- Trace file setup ---
  if (getenv("MALLOC_TRACE")) {
    if (!default_allocator) {
      default_allocator = a;
    }
    const char* file = getenv("MALLOC_TRACE");
    if (!file || std::strlen(file) == 0) {
      a->trace_out = stderr;
    } else {
      char filename[MAXPATHLEN];
      std::sprintf(filename, file, getpid());
      a->trace_out = std::fopen(filename, "w");
      std::setvbuf(a->trace_out, trace_buffer, _IOFBF, STATSIZE);
      if (!a->trace_out) {
        std::perror("fopen");
        std::fprintf(stderr, "cannot open trace file: %s, not tracing\n", file);
        a->trace_out = nullptr;
      }
    }
    if (a->trace_out && !a->stats_out) {
      a->stats_out = a->trace_out;
#if !defined(USE_LENNY_HACK)
      std::atexit(shutdown);
      atexit_added = true;
#endif
    }
  } else {
    a->trace_out = nullptr;
  }

  // --- Stats file setup ---
  a->statsfile    = nullptr;
  const char* statsfile = getenv("MALLOC_STATS");
  if (statsfile) {
    if (!default_allocator) {
      default_allocator = a;
    }
    if (std::strlen(statsfile) == 0) {
      a->stats_out = stderr;
    } else {
      a->statsfile = statsfile;
      a->stats_out = nullptr;
    }
    if ((a->stats_out || statsfile) && !atexit_added) {
#if !defined(USE_LENNY_HACK)
      std::atexit(shutdown);
#endif
    }
  } else {
    a->stats_out = nullptr;
  }

  a->dieing = false;
  return a;
#else
  return nullptr;
#endif // DISABLE_SCI_MALLOC
}

// ---------------------------------------------------------------------------
// alloc
// ---------------------------------------------------------------------------
[[nodiscard]] void*
Allocator::alloc(std::size_t size, const char* tag, int linenum)
{
  if (size > MEDIUM_THRESHOLD) {
    return alloc_big(size, tag, linenum);
  }
  if (size == 0) {
    return nullptr;
  }

  AllocBin* obj_bin = get_bin(size);

#ifndef DEBUG
  if (obj_bin->maxsize < size || size < obj_bin->minsize) {
    std::fprintf(stderr, "maxsize: %zu\n", obj_bin->maxsize);
    std::fprintf(stderr, "size: %zu\n",    size);
    AllocError("Bins messed up...");
  }
#endif

  lock();

  if (!obj_bin->free) {
    fill_bin(obj_bin);
  }

  Tag* obj      = obj_bin->free;
  obj_bin->free = obj->next;
  if (obj_bin->free) {
    obj_bin->free->prev = nullptr;
  }

  obj->hunk->ninuse++;
  obj->tag = tag;
#ifdef USE_TAG_LINENUM
  obj->linenum = linenum;
#endif
  obj->next = obj_bin->inuse;
  if (obj_bin->inuse) {
    obj_bin->inuse->prev = obj;
  }
  obj->prev      = nullptr;
  obj_bin->inuse = obj;
  obj->reqsize   = size;
  obj_bin->ninuse++;

  nalloc++;
  sizealloc += size;
  std::size_t bytes_inuse = (sizealloc >= sizefree) ? (sizealloc - sizefree) : 0;
  if (bytes_inuse > highwater_alloc) {
    highwater_alloc = bytes_inuse;
  }
  obj_bin->nalloc++;

  unlock();

  if (!lazy) {
    audit(obj, OBJFREE);
  }

  auto* data  = static_cast<char*>(static_cast<void*>(obj));
  data       += sizeof(Tag);
  auto* sent1  = reinterpret_cast<Sentinel*>(data);
  data        += sizeof(Sentinel);
  char* d      = data;
  data        += obj_maxsize(obj);
  auto* sent2  = reinterpret_cast<Sentinel*>(data);

  sent1->first_word = sent1->second_word =
  sent2->first_word = sent2->second_word = SENT_VAL_INUSE;

  if (strict) {
    unsigned int start = static_cast<unsigned int>(
      (obj->reqsize + sizeof(int)) / sizeof(int));
    for (auto* p = reinterpret_cast<unsigned int*>(d) + start;
         p < reinterpret_cast<unsigned int*>(sent2); p++) {
      *p++ = FILL_PATTERN;
    }
  }

  if (trace_out) {
#ifdef USE_TAG_LINENUM
    std::fprintf(trace_out, "A %p %zu (%s:%d)\n", static_cast<void*>(d), size, tag, linenum);
#else
    std::fprintf(trace_out, "A %p %zu (%s)\n",    static_cast<void*>(d), size, tag);
#endif
  }

  if (do_shutdown) {
    shutdown();
  }
  return static_cast<void*>(d);
}

// ---------------------------------------------------------------------------
// alloc_big
// ---------------------------------------------------------------------------
[[nodiscard]] void*
Allocator::alloc_big(std::size_t size, const char* tag,
                     [[maybe_unused]] int linenum)
{
  lock();

  Tag*        obj    = big_bin.free;
  std::size_t osize  = size + OVERHEAD;
  std::size_t maxsize = osize + (size >> 4);

  for (; obj != nullptr; obj = obj->next) {
    if (obj->hunk->len > osize && obj->hunk->len <= maxsize) {
      break;
    }
  }

  if (!obj) {
    // Optionally trim the free list if it's grown too large.
    int nfree = big_bin.ntotal - big_bin.ninuse;
    if (nfree >= 20) {
      Tag* cursor = big_bin.free;
      for (int i = 0; i < 10; i++) {
        cursor = cursor->next;
      }
      Tag* last    = cursor;
      Tag* to_free = cursor->next;
      while (to_free != nullptr) {
        Tag*    next = to_free->next;
        OSHunk* h    = to_free->hunk;
        nmunmap++;
        sizemunmap += h->len + sizeof(OSHunk);
        OSHunk::free(h);
        to_free = next;
        big_bin.ntotal--;
      }
      last->next = nullptr;
    }

    // Allocate a new hunk.
    std::size_t tsize  = sizeof(OSHunk) + OVERHEAD + size;
    std::size_t npages = (tsize + pagesize - 1) / pagesize;
    tsize              = npages * pagesize;
    tsize             -= sizeof(OSHunk);

    unsigned long offset = sizeof(OSHunk) % ALIGN;
    if (offset != 0) {
      offset = ALIGN - offset;
    }
    tsize -= offset;

    OSHunk* hunk = OSHunk::alloc(tsize, true, this);
    nmmap++;
    sizemmap += tsize + sizeof(OSHunk);
    std::size_t diffmmap = sizemmap - sizemunmap;
    if (diffmmap > highwater_mmap) {
      highwater_mmap = diffmmap;
    }

    obj      = static_cast<Tag*>(hunk->data);
    obj->bin = &big_bin;
    obj->tag = "never used (big object)";
#ifdef USE_TAG_LINENUM
    obj->linenum = 0;
#endif
    obj->next = big_bin.free;
    if (big_bin.free) {
      big_bin.free->prev = obj;
    }
    obj->prev    = nullptr;
    big_bin.free = obj;
    obj->hunk    = hunk;
    big_bin.ntotal++;

    // Set up sentinels in the newly created object.
    {
      auto* data  = static_cast<char*>(static_cast<void*>(obj));
      data       += sizeof(Tag);
      auto* s1    = reinterpret_cast<Sentinel*>(data);
      data       += sizeof(Sentinel);
      char* d     = data;
      data       += obj_maxsize(obj);
      auto* s2    = reinterpret_cast<Sentinel*>(data);

      s1->first_word = s1->second_word =
      s2->first_word = s2->second_word = SENT_VAL_FREE;

      if (strict) {
        for (auto* p = reinterpret_cast<unsigned int*>(d);
             p < reinterpret_cast<unsigned int*>(s2); p++) {
          *p++ = FILL_PATTERN;
        }
      }
    }
  }

  // Unlink from free list.
  if (obj->prev) {
    obj->prev->next = obj->next;
  } else {
    big_bin.free = obj->next;
  }
  if (obj->next) {
    obj->next->prev = obj->prev;
  }

  obj->hunk->ninuse++;
  obj->tag = tag;
#ifdef USE_TAG_LINENUM
  obj->linenum = linenum;
#endif
  obj->next = big_bin.inuse;
  obj->prev = nullptr;
  if (big_bin.inuse) {
    big_bin.inuse->prev = obj;
  }
  big_bin.inuse = obj;
  obj->reqsize  = size;
  big_bin.ninuse++;

  nalloc++;
  sizealloc += size;
  std::size_t bytes_inuse = sizealloc - sizefree;
  if (bytes_inuse > highwater_alloc) {
    highwater_alloc = bytes_inuse;
  }
  big_bin.nalloc++;

  unlock();

  if (!lazy) {
    audit(obj, OBJFREE);
  }

  auto* data  = static_cast<char*>(static_cast<void*>(obj));
  data       += sizeof(Tag);
  auto* sent1  = reinterpret_cast<Sentinel*>(data);
  data        += sizeof(Sentinel);
  char* d      = data;
  data        += obj_maxsize(obj);
  auto* sent2  = reinterpret_cast<Sentinel*>(data);

  sent1->first_word = sent1->second_word =
  sent2->first_word = sent2->second_word = SENT_VAL_INUSE;

  if (strict) {
    unsigned int start = static_cast<unsigned int>(
      (obj->reqsize + sizeof(int)) / sizeof(int));
    for (auto* p = reinterpret_cast<unsigned int*>(d) + start;
         p < reinterpret_cast<unsigned int*>(sent2); p++) {
      *p++ = FILL_PATTERN;
    }
  }

  if (trace_out) {
#ifdef USE_TAG_LINENUM
    std::fprintf(trace_out, "A %p %zu (%s:%d)\n", static_cast<void*>(d), size, tag, linenum);
#else
    std::fprintf(trace_out, "A %p %zu (%s)\n",    static_cast<void*>(d), size, tag);
#endif
  }

  if (do_shutdown) {
    shutdown();
  }
  return static_cast<void*>(d);
}

// ---------------------------------------------------------------------------
// realloc
// ---------------------------------------------------------------------------
[[nodiscard]] void*
Allocator::realloc(void* dobj, std::size_t newsize)
{
  if (!dobj) {
    return alloc(newsize, "realloc", 0);
  }

  auto* dd      = static_cast<char*>(dobj) - sizeof(Sentinel) - sizeof(Tag);
  auto* oldobj  = reinterpret_cast<Tag*>(dd);

  if (!lazy) {
    audit(oldobj, OBJFREEING);
  }

  if (oldobj->bin) {
    AllocBin* oldbin = get_bin(oldobj->bin->maxsize);
    if (newsize <= obj_maxsize(oldobj) && newsize >= oldbin->minsize) {
      size_t oldsize  = oldobj->reqsize;
      oldobj->reqsize = newsize;

      auto* data  = static_cast<char*>(static_cast<void*>(oldobj));
      data       += sizeof(Tag) + sizeof(Sentinel);
      char* d     = data;
      data       += obj_maxsize(oldobj);
      auto* sent2 = reinterpret_cast<Sentinel*>(data);

      if (strict) {
        unsigned int start = static_cast<unsigned int>(
          (oldobj->reqsize + sizeof(int)) / sizeof(int));
        for (auto* p = reinterpret_cast<unsigned int*>(d) + start;
             p < reinterpret_cast<unsigned int*>(sent2); p++) {
          *p++ = FILL_PATTERN;
        }
      }

      if (trace_out) {
#ifdef USE_TAG_LINENUM
        std::fprintf(trace_out, "R %p %zu %p %zu (%s:%d)\n",
                     dobj, oldsize, dobj, newsize, oldobj->tag, oldobj->linenum);
#else
        std::fprintf(trace_out, "R %p %zu %p %zu (%s)\n",
                     dobj, oldsize, dobj, newsize, oldobj->tag);
#endif
      }
      return dobj;
    }
  }

  void*       nobj    = alloc(newsize, "realloc", 0);
  std::size_t oldsize = oldobj->reqsize;
  std::size_t minsize = (newsize < oldsize) ? newsize : oldsize;
  std::memcpy(nobj, dobj, minsize);   // Note: memcpy(dest, src, n) — reversed vs bcopy
  free(dobj);

  if (trace_out) {
#ifdef USE_TAG_LINENUM
    std::fprintf(trace_out, "R %p %zu %p %zu (%s:%d)\n",
                 dobj, oldsize, nobj, newsize, oldobj->tag, oldobj->linenum);
#else
    std::fprintf(trace_out, "R %p %zu %p %zu (%s)\n",
                 dobj, oldsize, nobj, newsize, oldobj->tag);
#endif
  }

  return nobj;
}

// ---------------------------------------------------------------------------
// memalign
// ---------------------------------------------------------------------------
[[nodiscard]] void*
Allocator::memalign(std::size_t alignment, std::size_t size, const char* ctag)
{
  if (alignment <= 8) {
    return alloc(size, ctag, 0);
  }

  std::size_t asize = size + OVERHEAD + alignment;
  void*       orig_addr = alloc(asize, ctag, 0);
  auto*       addr  = static_cast<char*>(orig_addr);

  std::size_t misalign =
    (reinterpret_cast<std::size_t>(addr) + OVERHEAD) % alignment;
  misalign = (misalign == 0) ? 0 : (alignment - misalign);
  addr    += misalign;

  auto* tag     = reinterpret_cast<Tag*>(addr);
  addr         += sizeof(Tag);
  tag->bin      = nullptr;
  tag->next     = tag->prev = reinterpret_cast<Tag*>(orig_addr);
  tag->hunk     = nullptr;
  tag->reqsize  = size;
  tag->tag      = ctag;
#ifdef USE_TAG_LINENUM
  tag->linenum = 0;
#endif
  auto* sent1             = reinterpret_cast<Sentinel*>(addr);
  addr                   += sizeof(Sentinel);
  sent1->first_word       = sent1->second_word = SENT_VAL_INUSE;
  return static_cast<void*>(addr);
}

// ---------------------------------------------------------------------------
// free
// ---------------------------------------------------------------------------
void
Allocator::free(void* dobj) noexcept
{
  if (!dobj) {
    return;
  }

  auto* dd  = static_cast<char*>(dobj) - sizeof(Sentinel) - sizeof(Tag);
  auto* obj = reinterpret_cast<Tag*>(dd);

  // std::cout << "DEBUG: Allocator::free(" << dobj << ") tag=" << (obj->tag ? obj->tag : "null") << std::endl;

  if (!obj->bin) {
    // Allocated via memalign
    if (!lazy) {
      audit(obj, OBJMEMALIGNFREEING);
    }
    if (obj->next != obj->prev) {
      AllocError("Memalign tag inconsistency, or memory corrupt!\n");
    }
    free(static_cast<void*>(obj->prev));
    if (do_shutdown) {
      shutdown();
    }
    return;
  }

  if (trace_out) {
#ifdef USE_TAG_LINENUM
    std::fprintf(trace_out, "F %p %zu (%s:%d)\n",
                 dobj, obj->reqsize, obj->tag, obj->linenum);
#else
    std::fprintf(trace_out, "F %p %zu (%s)\n",
                 dobj, obj->reqsize, obj->tag);
#endif
  }

  if (!lazy) {
    audit(obj, OBJFREEING);
  }

  AllocBin* obj_bin = get_bin(obj->bin->maxsize);

  lock();
  nfree++;
  sizefree += obj->reqsize;
  obj_bin->nfree++;

  // Remove from inuse list.
  if (obj->next) { obj->next->prev = obj->prev; }
  if (obj->prev) {
    obj->prev->next = obj->next;
  } else {
    obj_bin->inuse = obj->next;
  }
  obj_bin->ninuse--;

  if (obj_bin == &big_bin && obj->reqsize > 50 * 1024 * 1024) {
    OSHunk* hunk = obj->hunk;
    nmunmap++;
    sizemunmap += hunk->len + sizeof(OSHunk);
    OSHunk::free(hunk);
    big_bin.ntotal--;
  } else {
    // Return to free list.
    obj->next = obj_bin->free;
    if (obj_bin->free) {
      obj_bin->free->prev = obj;
    }
    obj->prev     = nullptr;
    obj_bin->free = obj;

    auto* data  = static_cast<char*>(static_cast<void*>(obj));
    data       += sizeof(Tag);
    auto* sent1  = reinterpret_cast<Sentinel*>(data);
    data        += sizeof(Sentinel);
    char* d      = data;
    data        += obj_maxsize(obj);
    auto* sent2  = reinterpret_cast<Sentinel*>(data);

    sent1->first_word = sent1->second_word =
    sent2->first_word = sent2->second_word = SENT_VAL_FREE;

    if (strict) {
      for (auto* p = reinterpret_cast<unsigned int*>(d);
           p < reinterpret_cast<unsigned int*>(sent2); p++) {
        *p++ = FILL_PATTERN;
      }
    }
  }
  unlock();

  if (do_shutdown) {
    shutdown();
  }
}

// ---------------------------------------------------------------------------
// fill_bin
// ---------------------------------------------------------------------------
void
Allocator::fill_bin(AllocBin* bin)
{
  nfillbin++;
  if (bin->maxsize > MEDIUM_THRESHOLD) {
    AllocError("fill_bin not finished...");
    return;
  }

  std::size_t tsize  = (bin->maxsize + OVERHEAD + ALIGN - 1) & ~(ALIGN - 1);
  auto        nalloc = static_cast<unsigned int>(SMALLEST_ALLOCSIZE / tsize);
  if (nalloc < 1) {
    nalloc = 1;
  }
  std::size_t reqsize = nalloc * tsize;

  OSHunk* hunk = nullptr;
  void*   p    = nullptr;
  get_hunk(reqsize, hunk, p);

  for (unsigned int i = 0; i < nalloc; i++) {
    auto* t  = static_cast<Tag*>(p);
    t->bin   = bin;
    t->tag   = "never used";
#ifdef USE_TAG_LINENUM
    t->linenum = 0;
#endif
    t->next = bin->free;
    if (bin->free) {
      bin->free->prev = t;
    }
    t->prev   = nullptr;
    bin->free = t;
    t->hunk   = hunk;
    p         = static_cast<void*>(static_cast<char*>(p) + tsize);

    auto* data  = static_cast<char*>(static_cast<void*>(t));
    data       += sizeof(Tag);
    auto* sent1  = reinterpret_cast<Sentinel*>(data);
    data        += sizeof(Sentinel);
    char* d      = data;
    data        += t->bin->maxsize;
    auto* sent2  = reinterpret_cast<Sentinel*>(data);

    sent1->first_word = sent1->second_word =
    sent2->first_word = sent2->second_word = SENT_VAL_FREE;

    if (strict) {
      for (auto* fp = reinterpret_cast<unsigned int*>(d);
           fp < reinterpret_cast<unsigned int*>(sent2); fp++) {
        *fp++ = FILL_PATTERN;
      }
    }
  }
  bin->ntotal += nalloc;
}

// ---------------------------------------------------------------------------
// init_bin
// ---------------------------------------------------------------------------
void
Allocator::init_bin(AllocBin* bin,
                    std::size_t maxsize,
                    std::size_t minsize) noexcept
{
  bin->maxsize = maxsize;
  bin->minsize = minsize;
  bin->free    = nullptr;
  bin->inuse   = nullptr;
  bin->ninuse  = 0;
  bin->ntotal  = 0;
}

// ---------------------------------------------------------------------------
// audit helpers
// ---------------------------------------------------------------------------
static void
printObjectAllocMessage(Tag* obj)
{
#ifdef USE_TAG_LINENUM
  std::fprintf(stderr,
               "Object was allocated with this tag:\n%s at this line number:%d\n",
               obj->tag, obj->linenum);
#else
  std::fprintf(stderr,
               "Object was allocated with this tag:\n%s\n", obj->tag);
#endif
}

void
Allocator::audit(Tag* obj, int what)
{
  auto* data  = static_cast<char*>(static_cast<void*>(obj));
  data       += sizeof(Tag);
  auto* sent1  = reinterpret_cast<Sentinel*>(data);
  data        += sizeof(Sentinel);
  char* d      = data;

  if (what != OBJMEMALIGNFREEING) {
    data += obj_maxsize(obj);
  }
  auto* sent2 = reinterpret_cast<Sentinel*>(data);

  if (what == OBJFREE) {
    if (sent1->first_word != SENT_VAL_FREE || sent1->second_word != SENT_VAL_FREE) {
      if (sent1->first_word == SENT_VAL_INUSE) {
        printObjectAllocMessage(obj);
        AllocError("Object should be free, but is tagged as INUSE");
      } else {
        std::fprintf(stderr, "DEBUG: sent1: 0x%x 0x%x (expected 0x%x)\n", 
                     sent1->first_word, sent1->second_word, SENT_VAL_FREE);
        std::fprintf(stderr, "Free object has been corrupted within\n");
        std::fprintf(stderr, "the 8 bytes before the allocated region\n");
        printObjectAllocMessage(obj);
        AllocError("Freed object corrupt");
      }
    }
    if (sent2->first_word != SENT_VAL_FREE || sent2->second_word != SENT_VAL_FREE) {
      if (sent2->first_word == SENT_VAL_INUSE) {
        AllocError("Object should be free, but is tagged as INUSE (on tail only)");
      } else {
        std::fprintf(stderr, "Free object has been corrupted within\n");
        std::fprintf(stderr, "the 8 bytes following the allocated region\n");
        printObjectAllocMessage(obj);
        AllocError("Freed object corrupt");
      }
    }
  } else if (what == OBJFREEING || what == OBJINUSE || what == OBJMEMALIGNFREEING) {
    if (sent1->first_word != SENT_VAL_INUSE || sent1->second_word != SENT_VAL_INUSE) {
      if (sent1->first_word == SENT_VAL_FREE) {
        if (what == OBJFREEING) {
          std::fprintf(stderr, "Pointer (%p) was freed twice!\n", static_cast<void*>(d));
          printObjectAllocMessage(obj);
          AllocError("Freeing pointer twice");
        } else {
          printObjectAllocMessage(obj);
          AllocError("Object should be inuse, but is tagged as FREE");
        }
      } else {
        std::fprintf(stderr, "DEBUG: sent1: 0x%x 0x%x (expected 0x%x)\n", 
                     sent1->first_word, sent1->second_word, SENT_VAL_INUSE);
        std::fprintf(stderr, "Object has been corrupted within\n");
        std::fprintf(stderr, "the 8 bytes before the allocated region\n");
        printObjectAllocMessage(obj);
        AllocError("Memory Object corrupt");
      }
    }
    if (what != OBJMEMALIGNFREEING) {
      if (sent2->first_word != SENT_VAL_INUSE || sent2->second_word != SENT_VAL_INUSE) {
        if (sent2->first_word == SENT_VAL_FREE) {
          if (what == OBJFREEING) {
            std::fprintf(stderr, "Pointer (%p) was freed twice! (tail only?)\n",
                         static_cast<void*>(d));
            printObjectAllocMessage(obj);
            AllocError("Freeing pointer twice");
          } else {
            printObjectAllocMessage(obj);
            AllocError("Object should be inuse, but is tagged as FREE");
          }
        } else {
          std::fprintf(stderr, "Object has been corrupted within\n");
          std::fprintf(stderr, "the 8 bytes after the allocated region\n");
          printObjectAllocMessage(obj);
          AllocError("Memory Object corrupt");
        }
      }
    }
  }

  if (strict && (what == OBJFREEING || what == OBJINUSE)) {
    unsigned int start = static_cast<unsigned int>(
      (obj->reqsize + sizeof(int)) / sizeof(int));
    for (auto* p = reinterpret_cast<unsigned int*>(d) + start;
         p < reinterpret_cast<unsigned int*>(sent2); p++) {
      unsigned int p1 = *p++;
      if (p1 != FILL_PATTERN) {
        std::fprintf(stderr, "p1=0x%x (should be 0x%x)\n",
                     static_cast<int>(p1), static_cast<int>(FILL_PATTERN));
        std::fprintf(stderr, "Object has been corrupted immediately ");
        std::fprintf(stderr, "after the allocated region\n");
        printObjectAllocMessage(obj);
        AllocError("Memory Object corrupt");
      }
    }
  }

  if (strict && what == OBJFREE) {
    for (auto* p = reinterpret_cast<unsigned int*>(d);
         p < reinterpret_cast<unsigned int*>(sent2); p++) {
      unsigned int p1 = *p++;
      if (p1 != FILL_PATTERN) {
        std::fprintf(stderr, "Object has been written after free\n");
        printObjectAllocMessage(obj);
        AllocError("Write after free");
      }
    }
  }
}

// ---------------------------------------------------------------------------
// get_hunk
// ---------------------------------------------------------------------------
void
Allocator::get_hunk(std::size_t reqsize, OSHunk*& ret_hunk, void*& ret_p)
{
  OSHunk* hunk = nullptr;
  for (hunk = hunks; hunk != nullptr; hunk = hunk->next) {
    if (hunk->spaceleft >= reqsize) {
      break;
    }
  }

  if (!hunk) {
    std::size_t s = (reqsize > NORMAL_OS_ALLOC_SIZE) ? reqsize : NORMAL_OS_ALLOC_SIZE;
    hunk            = OSHunk::alloc(s, false, this);
    hunk->next      = hunks;
    hunks           = hunk;
    hunk->spaceleft = s;
    hunk->curr      = hunk->data;
    nmmap++;
    sizemmap += s + sizeof(OSHunk);
    std::size_t diffmmap = sizemmap - sizemunmap;
    if (diffmmap > highwater_mmap) {
      highwater_mmap = diffmmap;
    }
  }

  hunk->spaceleft -= reqsize;
  ret_p            = hunk->curr;
  hunk->curr       = static_cast<void*>(static_cast<char*>(hunk->curr) + reqsize);

  if (default_allocator && default_allocator->trace_out) {
    std::fprintf(default_allocator->trace_out,
                 "H %p %p %zu\n",
                 static_cast<void*>(hunk),
                 ret_p,
                 reqsize);
  }
  ret_hunk = hunk;
}

// ---------------------------------------------------------------------------
// PrintTag
// ---------------------------------------------------------------------------
void
PrintTag(void* dobj)
{
  auto* dd    = static_cast<char*>(dobj) - sizeof(Sentinel);
  auto* sent1 = reinterpret_cast<Sentinel*>(dd);
  dd         -= sizeof(Tag);
  auto* obj   = reinterpret_cast<Tag*>(dd);

#ifdef USE_TAG_LINENUM
  std::fprintf(stderr, "tag %p: allocated by: %s at %d\n",
               static_cast<void*>(obj), obj->tag, obj->linenum);
#else
  std::fprintf(stderr, "tag %p: allocated by: %s\n",
               static_cast<void*>(obj), obj->tag);
#endif
  std::fprintf(stderr, "requested object size: %zu bytes\n", obj->reqsize);
  std::fprintf(stderr, "maximum bin size: %zu bytes\n",      obj->bin->maxsize);
  std::fprintf(stderr, "range of object: %p - %zu\n",
               dobj,
               reinterpret_cast<std::size_t>(dobj) + obj->reqsize);
  std::fprintf(stderr,
               "range of object with overhead and sentinels: %p - %p\n",
               static_cast<void*>(obj),
               static_cast<void*>(reinterpret_cast<char*>(obj) + OVERHEAD));
  std::fprintf(stderr, "range of hunk: %zu - %zu\n",
               reinterpret_cast<std::size_t>(obj->hunk->data),
               reinterpret_cast<std::size_t>(obj->hunk->data) + obj->hunk->len);
  std::fprintf(stderr, "pre-sentinels: %x %x\n",
               sent1->first_word, sent1->second_word);

  if (sent1->first_word == SENT_VAL_FREE && sent1->second_word == SENT_VAL_FREE) {
    std::fprintf(stderr, "object should be free\n");
  } else if (sent1->first_word == SENT_VAL_INUSE && sent1->second_word == SENT_VAL_INUSE) {
    std::fprintf(stderr, "object should be inuse\n");
  } else {
    std::fprintf(stderr, "status of object is unknown - sentinels must be messed up\n");
  }
}

// ---------------------------------------------------------------------------
// GetGlobalStats (full accounting)
// ---------------------------------------------------------------------------
void
GetGlobalStats(Allocator*   a,
               std::size_t& nalloc,
               std::size_t& sizealloc,
               std::size_t& nfree,
               std::size_t& sizefree,
               std::size_t& nfillbin,
               std::size_t& nmmap,
               std::size_t& sizemmap,
               std::size_t& nmunmap,
               std::size_t& sizemunmap,
               std::size_t& highwater_alloc,
               std::size_t& highwater_mmap,
               std::size_t& bytes_overhead,
               std::size_t& bytes_free,
               std::size_t& bytes_fragmented,
               std::size_t& bytes_inuse,
               std::size_t& bytes_inhunks)
{
  if (!a) {
    nalloc = sizealloc = nfree = sizefree = nfillbin = 0;
    nmmap = sizemmap = nmunmap = sizemunmap          = 0;
    highwater_alloc = highwater_mmap                 = 0;
    bytes_overhead = bytes_free = bytes_fragmented = bytes_inuse = bytes_inhunks = 0;
    return;
  }

  a->lock();
  nalloc          = a->nalloc;
  sizealloc       = a->sizealloc;
  nfree           = a->nfree;
  sizefree        = a->sizefree;
  nfillbin        = a->nfillbin;
  nmmap           = a->nmmap;
  sizemmap        = a->sizemmap;
  nmunmap         = a->nmunmap;
  sizemunmap      = a->sizemunmap;
  highwater_alloc = a->highwater_alloc;
  highwater_mmap  = a->highwater_mmap;

  bytes_overhead = bytes_free = bytes_fragmented = bytes_inuse = bytes_inhunks = 0;

  for (int i = 0; i < static_cast<int>(NSMALL_BINS); i++) {
    account_bin(a, &a->small_bins[i], nullptr,
                bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);
  }
  for (int i = 0; i < static_cast<int>(NMEDIUM_BINS); i++) {
    account_bin(a, &a->medium_bins[i], nullptr,
                bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);
  }
  account_bin(a, &a->big_bin, nullptr,
              bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse);

  for (OSHunk* hunk = a->hunks; hunk != nullptr; hunk = hunk->next) {
    bytes_overhead += sizeof(OSHunk);
    bytes_inhunks  += hunk->spaceleft;
  }
  for (Tag* p = a->big_bin.free;  p != nullptr; p = p->next) { bytes_overhead += sizeof(OSHunk); }
  for (Tag* p = a->big_bin.inuse; p != nullptr; p = p->next) { bytes_overhead += sizeof(OSHunk); }
  bytes_overhead += a->mysize;

  a->unlock();
}

// ---------------------------------------------------------------------------
// GetGlobalStats (fast — no bin accounting)
// ---------------------------------------------------------------------------
void
GetGlobalStats(Allocator*   a,
               std::size_t& nalloc,
               std::size_t& sizealloc,
               std::size_t& nfree,
               std::size_t& sizefree,
               std::size_t& nfillbin,
               std::size_t& nmmap,
               std::size_t& sizemmap,
               std::size_t& nmunmap,
               std::size_t& sizemunmap,
               std::size_t& highwater_alloc,
               std::size_t& highwater_mmap)
{
  if (!a) {
    nalloc = sizealloc = nfree = sizefree = nfillbin = 0;
    nmmap = sizemmap = nmunmap = sizemunmap           = 0;
    highwater_alloc = highwater_mmap                  = 0;
    return;
  }

  a->lock();
  nalloc          = a->nalloc;
  sizealloc       = a->sizealloc;
  nfree           = a->nfree;
  sizefree        = a->sizefree;
  nfillbin        = a->nfillbin;
  nmmap           = a->nmmap;
  sizemmap        = a->sizemmap;
  nmunmap         = a->nmunmap;
  sizemunmap      = a->sizemunmap;
  highwater_alloc = a->highwater_alloc;
  highwater_mmap  = a->highwater_mmap;
  a->unlock();
}

// ---------------------------------------------------------------------------
// Bin stats
// ---------------------------------------------------------------------------
int
GetNbins(Allocator*)
{
  return static_cast<int>(NSMALL_BINS + NMEDIUM_BINS + 1);
}

void
GetBinStats(Allocator*   a,
            int          binno,
            std::size_t& minsize,
            std::size_t& maxsize,
            std::size_t& nalloc,
            std::size_t& nfree,
            std::size_t& ninlist)
{
  AllocBin* bin = nullptr;
  if (binno < static_cast<int>(NSMALL_BINS)) {
    bin = &a->small_bins[binno];
  } else if (binno < static_cast<int>(NSMALL_BINS + NMEDIUM_BINS)) {
    bin = &a->medium_bins[binno - static_cast<int>(NSMALL_BINS)];
  } else {
    bin = &a->big_bin;
  }

  a->lock();
  minsize = bin->minsize;
  maxsize = bin->maxsize;
  nalloc  = bin->nalloc;
  nfree   = bin->nfree;
  ninlist = bin->ntotal - bin->ninuse;
  a->unlock();
}

// ---------------------------------------------------------------------------
// DefaultAllocator
// ---------------------------------------------------------------------------
[[nodiscard]] Allocator*
DefaultAllocator()
{
  MakeDefaultAllocator();
  return default_allocator;
}

// ---------------------------------------------------------------------------
// Audit helpers
// ---------------------------------------------------------------------------
static void
audit_bin(Allocator* a, AllocBin* bin)
{
  for (Tag* p = bin->free; p != nullptr; p = p->next) {
    if (p->next && p->next->prev != p) {
      AllocError("Free list confused");
    }
    a->audit(p, OBJFREE);
  }
  for (Tag* p = bin->inuse; p != nullptr; p = p->next) {
    if (p->next && p->next->prev != p) {
      AllocError("Inuse list confused");
    }
    a->audit(p, OBJINUSE);
  }
}

void
AuditAllocator(Allocator* a)
{
  a->lock();
  for (int i = 0; i < static_cast<int>(NSMALL_BINS);  i++) { audit_bin(a, &a->small_bins[i]);  }
  for (int i = 0; i < static_cast<int>(NMEDIUM_BINS); i++) { audit_bin(a, &a->medium_bins[i]); }
  audit_bin(a, &a->big_bin);
  a->unlock();
}

void
AuditDefaultAllocator()
{
  AuditAllocator(default_allocator);
}

// ---------------------------------------------------------------------------
// DumpAllocator
// ---------------------------------------------------------------------------
static void
dump_bin(Allocator*, AllocBin* bin, FILE* fp)
{
  for (Tag* p = bin->inuse; p != nullptr; p = p->next) {
    auto* obj_ptr = static_cast<void*>(
      reinterpret_cast<char*>(p) + sizeof(Tag) + sizeof(Sentinel));
#ifdef USE_TAG_LINENUM
    std::fprintf(fp, "%p %zu %s:%d\n", obj_ptr, p->reqsize, p->tag, p->linenum);
#else
    std::fprintf(fp, "%p %zu %s\n",    obj_ptr, p->reqsize, p->tag);
#endif
  }
}

void
DumpAllocator(Allocator* a, const char* filename)
{
  if (!a) {
    std::printf("WARNING: In DumpAllocator: Allocator is nullptr.\n");
    std::printf("         Therefore no information to dump.");
    return;
  }

  FILE* fp = std::fopen(filename, "w");
  if (!fp) {
    std::perror("DumpAllocator fopen");
    std::exit(1);
  }

  std::fprintf(fp, "\n");
  a->lock();
  for (int i = 0; i < static_cast<int>(NSMALL_BINS);  i++) { dump_bin(a, &a->small_bins[i],  fp); }
  for (int i = 0; i < static_cast<int>(NMEDIUM_BINS); i++) { dump_bin(a, &a->medium_bins[i], fp); }
  dump_bin(a, &a->big_bin, fp);
  a->unlock();
  std::fclose(fp);
}

// ---------------------------------------------------------------------------

void
Allocator::noninline_unlock()
{
  unlock();
}

} // namespace Uintah

#endif // !defined(DISABLE_SCI_MALLOC)