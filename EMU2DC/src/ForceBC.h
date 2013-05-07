#ifndef __EMU2DC_FORCEBC_H__
#define __EMU2DC_FORCEBC_H__

#include <NodePArray.h>
#include <Types.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

namespace Emu2DC {
  
  class ForceBC {
  
  public:

    ForceBC();
    ~ForceBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes);

  private:

    void findBoundaryNodesInBox(const SCIRun::Vector& boxMin, 
                                const SCIRun::Vector& boxMax,
                                const NodePArray& nodes, 
                                NodePArray& boundaryNodes);

    void computeExtForceDensity(const SCIRun::Vector& extForce,
                                NodePArray& boundaryNodes);

    void initialize(std::string input_data_file);
    void computeExtForceDensity(const NodePArray& nodes,
	  	                const Array3& topLoc,
		                const Array3& botLoc,
                                const Array3& extForce);
    void findBoundaryNodes(const NodePArray& nodes,
	  	           const Array3& topLoc,
		           const Array3& botLoc,
		           NodePArray& topNodes,
			   NodePArray& botNodes);
    void sortNodes(const NodePArray& boundaryNodes,
		   NodePArray& sortedBoundaryNodes,
		   std::vector<double>& nodalSpan);

   // prevent copying
   ForceBC(const ForceBC& dyna);
   ForceBC& operator=(const ForceBC& dyna);

  }; // end class

} // end namespace
#endif

