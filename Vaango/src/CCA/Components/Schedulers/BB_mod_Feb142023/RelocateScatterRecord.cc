/*
 * The MIT License
 *
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

#include <CCA/Components/Schedulers/RelocateScatterRecord.h>

#include <iostream>

namespace Uintah {

std::ostream&
operator<<(std::ostream& out, const ScatterRecord& r)
{
  out.setf(std::ios::scientific, std::ios::floatfield);
  out.precision(4);
  out << " Scatter Record, matl: " << r.matl << " Level: " << r.level_index
      << " numParticles " << r.send_pset->numParticles()
      << " (Particle moving from Patch " << r.from_patch->getID()
      << ", to Patch " << r.to_patch->getID() << ")"
      << " vectorToNeighbor " << r.vector_to_neighbor;
  out.setf(std::ios::scientific, std::ios::floatfield);
  return out;
}

} // namespace Uintah
