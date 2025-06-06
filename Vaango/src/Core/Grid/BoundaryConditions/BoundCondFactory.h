/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_BOUND_COND_FACTORY_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_BOUND_COND_FACTORY_H__

#include <Core/Grid/BoundaryConditions/BoundCondBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

namespace Uintah {

class BoundCondFactory
{
public:
  // this function has a switch for all known BC_types
  static auto
  create(ProblemSpecP& ps, int& mat_id, const std::string face_label) -> BoundCondBaseSP;

  static auto
  customBC(int mat_id,
           const std::string face_label,
           double value,
           std::string label,
           std::string var) -> BoundCondBaseSP;

  static auto
  customBC(int mat_id,
           const std::string face_label,
           Vector value,
           std::string label,
           std::string var) -> BoundCondBaseSP;

  static auto
  customBC(int mat_id,
           const std::string face_label,
           std::string value,
           std::string label,
           std::string var) -> BoundCondBaseSP;
};

} // End namespace Uintah

#endif /* __CORE_GRID_BOUNDARYCONDITIONS_BOUND_COND_FACTORY_H__ */
