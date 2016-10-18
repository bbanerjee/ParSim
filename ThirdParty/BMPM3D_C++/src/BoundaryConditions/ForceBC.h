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

#ifndef __MATITI_FORCEBC_H__
#define __MATITI_FORCEBC_H__

#include <NodePArray.h>
#include <Types.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

namespace BrMPM {
  
  class ForceBC {
  
  public:

    ForceBC();
    ~ForceBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes);

  private:

    void findSurfaceNodesInBox(const SCIRun::Vector& boxMin, 
                                const SCIRun::Vector& boxMax,
                                const NodePArray& nodes, 
                                NodePArray& surfaceNodes);

    void findMaxVolume(NodePArray& surfaceNodes, double& maxVol);

    void computeExtForceDensity(const SCIRun::Vector& extForce,
                                NodePArray& surfaceNodes, double& maxvol);

    void initialize(std::string input_data_file);
    void computeExtForceDensity(const NodePArray& nodes,
	  	                const Point3D& topLoc,
		                const Point3D& botLoc,
                                const Vector3D& extForce);
    void findSurfaceNodes(const NodePArray& nodes,
	  	           const Point3D& topLoc,
		           const Point3D& botLoc,
		           NodePArray& topNodes,
			   NodePArray& botNodes);
    void sortNodes(const NodePArray& surfaceNodes,
		   NodePArray& sortedSurfaceNodes,
		   std::vector<double>& nodalSpan);

   // prevent copying
   ForceBC(const ForceBC& dyna);
   ForceBC& operator=(const ForceBC& dyna);

  }; // end class

} // end namespace
#endif

