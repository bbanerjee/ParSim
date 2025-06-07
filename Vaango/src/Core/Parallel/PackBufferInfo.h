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

#ifndef __CORE_PARALLEL_PACKBUFFERINFO_H__
#define __CORE_PARALLEL_PACKBUFFERINFO_H__

#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/BufferInfo.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Parallel/UintahMPI.h>
#include <Core/Util/RefCounted.h>

namespace Uintah {

class PackedBuffer : public RefCounted
{
public:
  PackedBuffer(int bytes)
    : d_buffer((void*)(scinew char[bytes]))
    , d_buffer_size(bytes)
  {
  }

  ~PackedBuffer()
  {
    delete[] (char*)d_buffer;
    d_buffer = nullptr;
  }

  void* getBuffer() { return d_buffer; }
  int getBufSize() { return d_buffer_size; }

private:
  void* d_buffer;
  int d_buffer_size;
};

class PackBufferInfo : public BufferInfo
{
public:
  PackBufferInfo();
  ~PackBufferInfo();

  PackBufferInfo(const PackBufferInfo&) = delete;
  PackBufferInfo(PackBufferInfo&&) = delete;

  PackBufferInfo& operator=(const PackBufferInfo&) = delete;
  PackBufferInfo& operator=(PackBufferInfo&&) = delete;

  void get_type(void*& output_buffer, int& output_count,
                MPI_Datatype& output_datatype, MPI_Comm comm);

  void get_type(void*& output_buffer, int& output_count,
                MPI_Datatype& output_datatype);

  void pack(MPI_Comm comm, int& out_count);

  void unpack(MPI_Comm comm, MPI_Status& status);

  // PackBufferInfo is to be an AfterCommuncationHandler object for the
  // MPI_CommunicationRecord template in MPIScheduler.cc.  After receive
  // requests have finished, then it needs to unpack what got received.
  void finishedCommunication(const ProcessorGroup* pg, MPI_Status& status)
  {
    unpack(pg->getComm(), status);
  }

private:
  PackedBuffer* d_packedBuffer{ nullptr };
};
} // namespace Uintah

#endif //__CORE_PARALLEL_PACKBUFFERINFO_H__
