/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Util/DebugStream.h>

namespace Uintah {

static DebugStream dbg("GeometryPiece", false);

void
GeometryPiece::outputProblemSpec(ProblemSpecP& ps) const {
  ProblemSpecP child_ps = ps->appendChild(getType().c_str());

  if (d_nameSet) {
    child_ps->setAttribute("label", d_name);

    if (d_firstOutput) {
      // If geom obj is named, then only output data the first time.
      dbg << "GP::outputProblemSpec(): Full description of: " << d_name
          << " -- " << getType() << "\n";
      outputHelper(child_ps);
      d_firstOutput = false;

    } else {
      dbg << "GP::outputProblemSpec(): Reference to: " << d_name << " -- "
          << getType() << "\n";
    }

  } else {
    dbg << "GP::outputProblemSpec(): Full Description Of: " << d_name << " -- "
        << getType() << "\n";
    // If no name, then always print out all data.
    outputHelper(child_ps);
  }
}

}  // end namespace Uintah