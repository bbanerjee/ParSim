/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2014-2023 Parresia Research Limited
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

#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/BufferInfo.h>
#include <Core/Util/Assert.h>
#include <Core/Util/RefCounted.h>

namespace Uintah {

BufferInfo::~BufferInfo() noexcept(false)
{
  if (d_free_datatype) {
    ASSERT(d_datatype != MPI_DATATYPE_NULL);
    ASSERT(d_datatype != MPI_INT);
    ASSERT(d_datatype != MPI_DOUBLE);
    Uintah::MPI::Type_free(&d_datatype);
    d_datatype = MPI_DATATYPE_NULL;
  }

  for (unsigned int i = 0; i < d_datatypes.size(); i++) {
    if (d_free_datatypes[i]) {
      ASSERT(d_datatypes[i] != MPI_DATATYPE_NULL);
      ASSERT(d_datatypes[i] != MPI_INT);
      ASSERT(d_datatypes[i] != MPI_DOUBLE);
      Uintah::MPI::Type_free(&d_datatypes[i]);
      d_datatypes[i] = MPI_DATATYPE_NULL;
    }
  }

  if (d_sendlist) {
    delete d_sendlist;
    d_sendlist = nullptr;
  }
}

unsigned int
BufferInfo::count() const
{
  return d_datatypes.size();
}

void
BufferInfo::add(void* startbuf, int count, MPI_Datatype datatype,
                bool free_datatype)
{
  ASSERT(!d_have_datatype);
  d_startbufs.push_back(startbuf);
  d_counts.push_back(count);
  d_datatypes.push_back(datatype);
  d_free_datatypes.push_back(free_datatype);
}

void
BufferInfo::get_type(void*& out_buf, int& out_count, MPI_Datatype& out_datatype)
{
  ASSERT(count() > 0);
  if (!d_have_datatype) {
    if (count() == 1) {
      d_buffer = d_startbufs[0];
      d_count = d_counts[0];
      d_datatype = d_datatypes[0];
      d_free_datatype = false; // Will get freed with array
    } else {
      //MPI_Type_create_struct allows for multiple things to be sent in a single message.
      std::vector<MPI_Aint> offsets(count());
      offsets[0] = 0;
      for (unsigned int i = 0; i < d_startbufs.size(); i++) {
        // Find how far this address is displaced from the first address.
        // From MPI's point of view, it won't know whether these offsets are from the same array or
        // from entirely different variables.  We'll take advantage of that second point.
        // It also appears to allow for negative displacements.
        offsets[i] =
          (MPI_Aint)((char*)d_startbufs.at(i) - (char*)d_startbufs.at(0));
      }
      Uintah::MPI::Type_create_struct(count(), &d_counts[0], &offsets[0], 
                                      &d_datatypes[0], &d_datatype);
      Uintah::MPI::Type_commit(&d_datatype);
      d_buffer = d_startbufs[0];
      d_count = 1;
      d_free_datatype = true;
    }
    d_have_datatype = true;
  }
  out_buf = d_buffer;
  out_count = d_count;
  out_datatype = d_datatype;
}

Sendlist::~Sendlist()
{
  if (d_obj && d_obj->removeReference()) {
    delete d_obj;
    d_obj = nullptr;
  }

  // Avoid recursive deletion
  Sendlist* p = d_next;
  while (p) {
    if (p->d_obj->removeReference())
      delete p->d_obj;
    Sendlist* n = p->d_next;
    p->d_next = 0; // So that DTOR won't recurse...
    p->d_obj = 0;
    delete p;
    p = n;
  }
}

void
BufferInfo::addSendlist(RefCounted* obj)
{
  obj->addReference();
  d_sendlist = scinew Sendlist(d_sendlist, obj);
}

Sendlist*
BufferInfo::takeSendlist()
{
  Sendlist* rtn = d_sendlist;
  d_sendlist = nullptr; // They are now responsible for freeing...
  return rtn;
}

} // end namespace Uintah