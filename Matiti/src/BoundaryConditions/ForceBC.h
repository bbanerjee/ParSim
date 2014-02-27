#ifndef __MATITI_FORCEBC_H__
#define __MATITI_FORCEBC_H__

#include <BoundaryConditions/LoadBC.h>

namespace Matiti {
  
  class ForceBC : public LoadBC {
  
  public:

    ForceBC();
    ~ForceBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elems);

  protected:

    void computeExtForceDensity(const SCIRun::Vector& extForce,
                                NodePArray& surfaceNodes, 
                                double& maxvol);

  private:

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

