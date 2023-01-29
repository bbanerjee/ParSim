/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <Core/Exceptions/InternalError.h>
#include <Core/Parallel/PackBufferInfo.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/Assert.h>
#include <Core/Util/RefCounted.h>

#include <iostream>
#include <sstream>
#include <string.h>
#include <zlib.h>

namespace Uintah {

PackBufferInfo::PackBufferInfo()
  : BufferInfo()
{
}

PackBufferInfo::~PackBufferInfo()
{
  if (d_packedBuffer && d_packedBuffer->removeReference()) {
    delete d_packedBuffer;
    d_packedBuffer = nullptr;
  }
}

void
PackBufferInfo::get_type(void*& out_buf, int& out_count,
                         MPI_Datatype& out_datatype, MPI_Comm comm)
{
  ASSERT(count() > 0);
  if (!d_have_datatype) {
    int packed_size;
    int total_packed_size = 0;
    for (unsigned int i = 0; i < d_startbufs.size(); i++) {
      if (d_counts[i] > 0) {
        MPI_Pack_size(d_counts[i], d_datatypes[i], comm, &packed_size);
        total_packed_size += packed_size;
      }
    }

    d_packedBuffer = scinew PackedBuffer(total_packed_size);
    d_packedBuffer->addReference();

    d_datatype = MPI_PACKED;
    d_count = total_packed_size;
    d_buffer = d_packedBuffer->getBuffer();
    d_have_datatype = true;
  }
  out_buf = d_buffer;
  out_count = d_count;
  out_datatype = d_datatype;
}

void
PackBufferInfo::get_type(void*&, int&, MPI_Datatype&)
{
  // Should use other overload for a PackBufferInfo
  std::ostringstream out;
  out << "get_type(void*&, int&, MPI_Datatype&) should not be "
      << "called on PackBufferInfo objects";
  throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
}

void
PackBufferInfo::pack(MPI_Comm comm, int& out_count)
{
  ASSERT(d_have_datatype);

  int position = 0;
  int bufsize = d_packedBuffer->getBufSize();

  unsigned int bufIndex = 0;
  for (auto& startbuf : d_startbufs) {
    unsigned int count = d_counts[bufIndex];
    if (count > 0) {
      // pack into a contigious buffer
      MPI_Pack(startbuf, count, d_datatypes[bufIndex], d_buffer, bufsize,
               &position, comm);
    }
    bufIndex++;
  }

  out_count = position;

  // When it is all packed, only the buffer necessarily needs to be kept
  // around until after it is sent.
  delete d_sendlist;
  d_sendlist = nullptr;
  addSendlist(d_packedBuffer);
}

void
PackBufferInfo::unpack(MPI_Comm comm, [[maybe_unused]] MPI_Status& status)
{
  ASSERT(d_have_datatype);

  unsigned long bufsize = d_packedBuffer->getBufSize();

  int position = 0;
  unsigned int bufIndex = 0;
  for (auto& startbuf : d_startbufs) {
    unsigned int count = d_counts[bufIndex];
    if (count > 0) {
      MPI_Unpack(d_buffer, bufsize, &position, startbuf, count,
                 d_datatypes[bufIndex], comm);
    }
  }
}

} // end namespace Uintah
