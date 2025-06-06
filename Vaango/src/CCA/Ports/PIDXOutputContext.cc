/*
 * The MIT License
 *
 * Copyright (c) 2013-2015 The University of Utah
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

#include <CCA/Ports/PIDXOutputContext.h>

#if HAVE_PIDX

namespace Uintah {

PIDXOutputContext::PIDXOutputContext() {}

PIDXOutputContext::~PIDXOutputContext()
{
  PIDX_close_access(this->access);
}

void
PIDXOutputContext::initialize(std::string filename,
                              unsigned int timeStep,
                              int globalExtents[3],
                              MPI_Comm comm)
{

  this->filename = filename;
  this->timestep = timeStep;
  this->comm     = comm;

  PIDX_point global_bounding_box;
  PIDX_create_access(&(this->access));
  PIDX_set_mpi_access(this->access, this->comm);

  PIDX_set_point_5D(global_bounding_box,
                    globalExtents[0],
                    globalExtents[1],
                    globalExtents[2],
                    1,
                    1);

  PIDX_file_create(filename.c_str(), PIDX_MODE_CREATE, access, &(this->file));

  int64_t restructured_box_size[5] = { 64, 64, 64, 1, 1 };
  PIDX_set_restructuring_box(file, restructured_box_size);

  // PIDX_set_resolution(this->file, 0, 2);
  PIDX_set_dims(this->file, global_bounding_box);
  PIDX_set_current_time_step(this->file, timeStep);
  PIDX_set_block_size(this->file, 16);
  PIDX_set_block_count(this->file, 128);

  PIDX_set_compression_type(this->file, PIDX_CHUNKING_ZFP);
  PIDX_set_lossy_compression_bit_rate(this->file, 8);
}

} // end namespace Uintah

#endif // HAVE_PIDX
