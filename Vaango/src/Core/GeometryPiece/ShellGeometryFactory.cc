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

#include <Core/GeometryPiece/CylinderShellPiece.h>
#include <Core/GeometryPiece/PlaneShellPiece.h>
#include <Core/GeometryPiece/ShellGeometryFactory.h>
#include <Core/GeometryPiece/SphereShellPiece.h>
// #include <Core/GeometryPiece/GUVSphereShellPiece.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <iostream>
#include <string>

namespace Uintah {

GeometryPieceP
ShellGeometryFactory::create(ProblemSpecP& ps) {
  std::string go_type = ps->getNodeName();

  if (go_type == PlaneShellPiece::TYPE_NAME) {
    return std::make_shared<PlaneShellPiece>(ps);
  } else if (go_type == SphereShellPiece::TYPE_NAME) {
    return std::make_shared<SphereShellPiece>(ps);
  } else if (go_type == CylinderShellPiece::TYPE_NAME) {
    return std::make_shared<CylinderShellPiece>(ps);
  }
  //  else if (go_type == GUVSphereShellPiece::TYPE_NAME ) {
  //    return std::make_unique<GUVSphereShellPiece>(ps);
  //  }
  return nullptr;
}

}  // end namespace Uintah
