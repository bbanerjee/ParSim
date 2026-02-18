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

#include <sci_defs/malloc_defs.h>

#include <Core/Malloc/mem_init.h>

#include <Core/Malloc/Allocator.h>

#ifdef MALLOC_TRACE
#include "MallocTraceOff.h"
#endif

#include <Core/Malloc/AllocPriv.h>
#include <new>
#include <cstdint>
#include <cassert>

using namespace Uintah;

#ifdef DISABLE_SCI_MALLOC

// These stubs are needed when your code uses these functions but
// DISABLE_SCI_MALLOC is set.
namespace Uintah {
const char*
AllocatorSetDefaultTagNew(const char* /*tag*/)
{
  return "AllocatorSetDefaultTagNew::NOT IMPLEMENTED.  DISABLE_SCI_MALLOC is "
         "set";
}

void
AllocatorResetDefaultTagNew()
{
}

const char*
AllocatorSetDefaultTag(const char* /*tag*/)
{
  return "AllocatorSetDefaultTag::NOT IMPLEMENTED.  DISABLE_SCI_MALLOC is set";
}

void
AllocatorResetDefaultTag()
{
}
int
AllocatorSetDefaultTagLineNumber(int line_number)
{
  return line_number;
}
void
AllocatorResetDefaultTagLineNumber()
{
}

} // namespace Uintah
#ifndef MALLOC_TRACE
void*
operator new(size_t size, Allocator*, const char*, int)
{
  void* mem = new char[size];
#ifdef INITIALIZE_MEMORY
  for (unsigned int i = 0; i < size; i++) {
    static_cast<unsigned char*>(mem)[i] = MEMORY_INIT_NUMBER;
  }
#endif
  return mem;
}

void*
operator new[](size_t size, Allocator*, const char*, int)
{
  void* mem = new char[size];

#ifdef INITIALIZE_MEMORY
  for (unsigned int i = 0; i < size; i++) {
    static_cast<unsigned char*>(mem)[i] = MEMORY_INIT_NUMBER;
  }
#endif
  return mem;
}

void
operator delete(void* ptr, Allocator*, const char*, int)
{
  free(ptr);
}

void
operator delete[](void* ptr, Allocator*, const char*, int)
{
  free(ptr);
}

#endif
#else // ifdef DISABLE_SCI_MALLOC = 0

static const char* default_new_tag       = "Unknown - operator new";
static const char* default_new_array_tag = "Unknown - operator new[]";

// the line number us an optional tag (on if configured with
// --enable-scinew-line-numbers) that can also show some information (like an
// interation or timestep) for each tag
int default_tag_line_number = 0;

namespace Uintah {
const char*
AllocatorSetDefaultTagNew(const char* tag)
{
  const char* old       = default_new_tag;
  default_new_tag       = tag;
  default_new_array_tag = tag;
  return old;
}

void
AllocatorResetDefaultTagNew()
{
  default_new_tag       = "Unknown - operator new";
  default_new_array_tag = "Unknown - operator new[]";
}

const char*
AllocatorSetDefaultTag(const char* tag)
{
  AllocatorSetDefaultTagMalloc(tag);
  return AllocatorSetDefaultTagNew(tag);
}

void
AllocatorResetDefaultTag()
{
  AllocatorResetDefaultTagMalloc();
  AllocatorResetDefaultTagNew();
}
int
AllocatorSetDefaultTagLineNumber(int line_number)
{
  int old_num             = default_tag_line_number;
  default_tag_line_number = line_number;
  return old_num;
}
void
AllocatorResetDefaultTagLineNumber()
{
  default_tag_line_number = 0;
}
} // namespace Uintah

// Helper to resolve allocator, avoiding repeated code
[[nodiscard]] static inline Allocator*
resolve_allocator(Allocator* a)
{
  if (!a) {
    if (!default_allocator) {
      MakeDefaultAllocator();
    }
    a = default_allocator;
  }
  return a;
}

// Helper to initialize memory (only active when INITIALIZE_MEMORY is defined)
static inline void
initialize_memory([[maybe_unused]] void* mem,
                  [[maybe_unused]] std::size_t size) noexcept
{
#ifdef INITIALIZE_MEMORY
  std::fill_n(static_cast<unsigned char*>(mem), size, MEMORY_INIT_NUMBER);
#endif
}

[[nodiscard]] static inline Allocator*
default_alloc()
{
  if (!default_allocator) {
    MakeDefaultAllocator();
  }
  return default_allocator;
}

// --- Standard new/delete ---

void*
operator new(std::size_t size)
{
  void* mem =
    default_alloc()->alloc(size, default_new_tag, default_tag_line_number);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  //printf("alignment: %zu\n", reinterpret_cast<std::uintptr_t>(mem) % 32);
  return mem;
}

void*
operator new[](std::size_t size)
{
  void* mem = default_alloc()->alloc(
    size, default_new_array_tag, default_tag_line_number);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  return mem;
}

// --- Nothrow variants ---

void*
operator new(std::size_t size, const std::nothrow_t&) noexcept
{
  void* mem = default_alloc()->alloc(
    size, "unknown - nothrow operator new", default_tag_line_number);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  return mem;
}

void*
operator new[](std::size_t size, const std::nothrow_t&) noexcept
{
  void* mem = default_alloc()->alloc(
    size, "unknown - nothrow operator new[]", default_tag_line_number);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  return mem;
}

// --- Aligned new/delete (C++17) ---

#if __cpp_aligned_new >= 201606L

void*
operator new(std::size_t size, std::align_val_t al)
{
  void* mem = default_alloc()->memalign(static_cast<std::size_t>(al), size, default_new_tag);
  initialize_memory(mem, size);
  return mem;
}

void*
operator new[](std::size_t size, std::align_val_t al)
{
  void* mem = default_alloc()->memalign(static_cast<std::size_t>(al), size, default_new_array_tag);
  initialize_memory(mem, size);
  return mem;
}

void*
operator new(std::size_t size, std::align_val_t al, const std::nothrow_t&) noexcept
{
  void* mem = default_alloc()->memalign(static_cast<std::size_t>(al), size, "unknown - nothrow aligned operator new");
  initialize_memory(mem, size);
  return mem;
}

void*
operator new[](std::size_t size, std::align_val_t al, const std::nothrow_t&) noexcept
{
  void* mem = default_alloc()->memalign(static_cast<std::size_t>(al), size, "unknown - nothrow aligned operator new[]");
  initialize_memory(mem, size);
  return mem;
}

void
operator delete(void* ptr, std::align_val_t) noexcept
{
  if (ptr) {
    default_alloc()->free(ptr);
  }
}

void
operator delete[](void* ptr, std::align_val_t) noexcept
{
  if (ptr) {
    default_alloc()->free(ptr);
  }
}

void
operator delete(void* ptr, std::size_t size, std::align_val_t al) noexcept
{
  if (ptr) {
    operator delete(ptr, al);
  }
  (void)size;
}

void
operator delete[](void* ptr, std::size_t size, std::align_val_t al) noexcept
{
  if (ptr) {
    operator delete[](ptr, al);
  }
  (void)size;
}

#endif

// --- Standard delete ---

void
operator delete(void* ptr) noexcept
{
  if (ptr) {
    default_alloc()->free(ptr);
  }
}

void
operator delete[](void* ptr) noexcept
{
  if (ptr) {
    default_alloc()->free(ptr);
  }
}

// --- Allocator-tagged new/delete ---

void*
operator new(std::size_t size, Allocator* a, const char* tag, int linenum)
{
  a         = resolve_allocator(a);
  void* mem = a->alloc(size, tag, linenum);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  return mem;
}

void*
operator new[](std::size_t size, Allocator* a, const char* tag, int linenum)
{
  a         = resolve_allocator(a);
  void* mem = a->alloc(size, tag, linenum);
  assert(reinterpret_cast<std::uintptr_t>(mem) % 32 == 0);
  initialize_memory(mem, size);
  return mem;
}

void
operator delete(void* ptr,
                Allocator* a,
                [[maybe_unused]] const char* tag,
                [[maybe_unused]] int linenum) noexcept
{
  if (ptr) {
    resolve_allocator(a)->free(ptr);
  }
}

void
operator delete[](void* ptr,
                  Allocator* a,
                  [[maybe_unused]] const char* tag,
                  [[maybe_unused]] int linenum) noexcept
{
  if (ptr) {
    resolve_allocator(a)->free(ptr);
  }
}

// Fix for -Wsized-deallocation warnings: sized deallocation overloads
void
operator delete(void* ptr, std::size_t size) noexcept
{
  if (ptr) {
    operator delete(ptr); // delegate to your existing unsized overload
  }
  (void)size;
}

void
operator delete[](void* ptr, std::size_t size) noexcept
{
  if (ptr) {
    operator delete[](ptr); // delegate to your existing unsized overload
  }
  (void)size;
}

#ifdef MALLOC_TRACE
#include "MallocTraceOn.h"
#endif

#endif // ifdef DISABLE_SCI_MALLOC
