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

#pragma once

#include <sci_defs/malloc_defs.h>
#include <sci_defs/thread_defs.h>

#include <cstddef>
#include <cstdio>

#ifdef SCI_PTHREAD
#include <pthread.h>
#else
#if !defined(SCI_NOTHREAD) && !defined(_WIN32)
#error "No lock implementation for this architecture"
#endif
#endif

namespace Uintah {

struct OSHunk;

// ---------------------------------------------------------------------------
// Sentinel — guards either side of every allocation to detect corruption.
// ---------------------------------------------------------------------------
struct Sentinel
{
  unsigned int first_word;
  unsigned int second_word;
};

struct AllocBin;

// ---------------------------------------------------------------------------
// Tag — header placed immediately before each allocation's sentinel pair.
// ---------------------------------------------------------------------------
struct Tag
{
  AllocBin*   bin;
  const char* tag;
#ifdef USE_TAG_LINENUM
  int         linenum;
  int         pad;     // Pad to make sizeof(Tag) + sizeof(Sentinel) a multiple of 32
#else
  long long   pad;     // Pad to make sizeof(Tag) + sizeof(Sentinel) a multiple of 32
#endif
  Tag*        next;
  Tag*        prev;
  OSHunk*     hunk;
  std::size_t reqsize;
};

// ---------------------------------------------------------------------------
// AllocBin — a size-class bucket holding free and in-use object lists.
// ---------------------------------------------------------------------------
struct AllocBin
{
  Tag*        free;
  Tag*        inuse;
  std::size_t maxsize;
  std::size_t minsize;
  int         ninuse;
  int         ntotal;
  std::size_t nalloc;
  std::size_t nfree;
};

// ---------------------------------------------------------------------------
// Allocator
// ---------------------------------------------------------------------------
struct Allocator
{
  // --- Locking ------------------------------------------------------------
#ifdef SCI_PTHREAD
  pthread_mutex_t the_lock;
#endif

  void        initlock() noexcept;
  inline void lock()    noexcept;
  inline void unlock()  noexcept;
  void        noninline_unlock();

#ifdef SCI_PTHREAD
  inline void rlock() noexcept;

  // These members deal with recursive locking behaviour needed to work
  // around glibc mutex/fork interaction bugs (see implementation for
  // details).  Initialised in initlock().
  bool      use_rlock;
  pthread_t owner;
  bool      owner_initialized;
  int       lock_count;
#endif

  // --- Core allocation / deallocation -------------------------------------
  [[nodiscard]] void* alloc_big(std::size_t size,
                                const char* tag,
                                int         linenum);

  [[nodiscard]] void* memalign(std::size_t alignment,
                               std::size_t size,
                               const char* tag);

  [[nodiscard]] void* alloc(std::size_t size,
                            const char* tag,
                            int         linenum);

#ifdef MALLOC_TRACE
#include <MallocTraceOff.h>
#endif

  void  free(void* ptr) noexcept;

  [[nodiscard]] void* realloc(void* ptr, std::size_t size);

#ifdef MALLOC_TRACE
#include <MallocTraceOn.h>
#endif

  // --- Configuration ------------------------------------------------------
  int   strict;
  int   lazy;
  FILE* trace_out;
  FILE* stats_out;
  const char* statsfile;    // changed from char* — never modified after set
  OSHunk* hunks;

  // --- Bin tables ---------------------------------------------------------
  AllocBin* small_bins;
  AllocBin* medium_bins;
  AllocBin  big_bin;

  // --- Internal helpers ---------------------------------------------------
  [[nodiscard]] inline AllocBin* get_bin(std::size_t size) noexcept;

  void fill_bin(AllocBin* bin);

  void get_hunk(std::size_t reqsize, OSHunk*& ret_hunk, void*& ret_p);

  void init_bin(AllocBin*   bin,
                std::size_t maxsize,
                std::size_t minsize) noexcept;

  void audit(Tag* obj, int what);

  [[nodiscard]] inline std::size_t obj_maxsize(Tag* t) noexcept;

  // --- Statistics ---------------------------------------------------------
  std::size_t nalloc;
  std::size_t nfree;
  std::size_t sizealloc;
  std::size_t sizefree;

  std::size_t nfillbin;
  std::size_t nmmap;
  std::size_t sizemmap;
  std::size_t nmunmap;
  std::size_t sizemunmap;

  std::size_t highwater_alloc;
  std::size_t highwater_mmap;

  std::size_t mysize;
  std::size_t pagesize;

  bool dieing;
};

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------
[[noreturn]] void AllocError(const char* msg);

} // namespace Uintah