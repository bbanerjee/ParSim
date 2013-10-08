#ifndef __MATITI_FORCEBC_H__
#define __MATITI_FORCEBC_H__

#include <NodePArray.h>
#include <Types.h>
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

namespace Matiti {
  
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

