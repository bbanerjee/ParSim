/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef __CCA_COMPONENTS_SCHEDULERS_SENDSTATE_H__
#define __CCA_COMPONENTS_SCHEDULERS_SENDSTATE_H__

#include <Core/Grid/Ghost.h>
#include <Core/Grid/Variables/PSPatchMatlGhostRange.h>

#include <map>

namespace Uintah {

class Patch;
class ParticleSubset;

class SendState
{
public:
  SendState() = default;
  ~SendState();

  // disable copy, assignment, and move
  SendState(const SendState&) = delete;
  SendState&
  operator=(const SendState&) = delete;
  SendState(SendState&&)      = delete;
  SendState&
  operator=(SendState&&) = delete;

  ParticleSubset*
  find_sendset(int dest,
               const Patch*,
               int matl,
               IntVector low,
               IntVector high,
               int dwid = 0) const;
  void
  add_sendset(ParticleSubset* pset,
              int dest,
              const Patch*,
              int matl,
              IntVector low,
              IntVector high,
              int dwid = 0);

  void
  reset(); // Clears out all sendsets...

  void
  print();

private:
  using map_type = std::map<std::pair<PSPatchMatlGhostRange, int>, ParticleSubset*>;
  map_type d_send_subsets;
};
} // End namespace Uintah

#endif //__CCA_COMPONENTS_SCHEDULERS_SENDSTATE_H__
