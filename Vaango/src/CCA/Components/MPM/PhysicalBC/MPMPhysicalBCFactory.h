/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#ifndef __Uintah_MPM_MPMPhysicalBCFactory__
#define __Uintah_MPM_MPMPhysicalBCFactory__

#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/Grid.h>
#include <CCA/Components/MPM/MPMFlags.h>

#include <memory>
#include <vector>


namespace Uintah {

  using MPMPhysicalBC_SP = std::shared_ptr<MPMPhysicalBC>;

  class MPMPhysicalBCFactory
  {
  public:
    static void create(const ProblemSpecP& ps, const GridP& grid, const MPMFlags* flags);
    static void clean(); // delete all mpmPhysicalBCs
    static std::vector<MPMPhysicalBC_SP> mpmPhysicalBCs;
  };
} // End namespace Uintah


#endif /* __Uintah_MPM_MPMPhysicalBCFactory__ */

