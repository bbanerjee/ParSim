/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MATITI_DISPLACEMENT_BC_H__
#define __MATITI_DISPLACEMENT_BC_H__

#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

namespace Matiti {
  
  class DisplacementBC {
  
  public:

    DisplacementBC();
    ~DisplacementBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes);
    void applyDisplacementBC();

  private:

    void initializeDispBCSurfaceNodes(const SCIRun::Vector& boxMin, 
                                      const SCIRun::Vector& boxMax,
                                      NodePArray& nodes);
  private:

    NodePArray d_surface_nodes;   // List of nodes to which these BCs are to be applied
    bool d_x_flag;     // Apply bc in the global x-direction if true
    bool d_y_flag;     // Apply bc in the global y-direction if true
    bool d_z_flag;     // Apply bc in the global z-direction if true
    double d_x_value;  // Value of bc in x-direction 
    double d_y_value;  // Value of bc in y-direction
    double d_z_value;  // Value of bc in z-direction

    // prevent copying
    DisplacementBC(const DisplacementBC& bc);
    DisplacementBC& operator=(const DisplacementBC& bc);

  }; // end class

} // end namespace
#endif

