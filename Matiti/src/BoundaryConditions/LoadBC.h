#ifndef __MATITI_LOAD_BC_H__
#define __MATITI_LOAD_BC_H__

#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>
#include <Types/Types.h>
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

namespace Matiti {
  
  class LoadBC {
  
  public:

    LoadBC();
    ~LoadBC();

    virtual void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elems) = 0;

  protected:

    void findSurfaceNodesInBox(const SCIRun::Vector& boxMin, 
                               const SCIRun::Vector& boxMax,
                               const NodePArray& nodes, 
                               NodePArray& surfaceNodes);

    void findMaxVolume(NodePArray& surfaceNodes, double& maxVol);

  private:

   // prevent copying
   LoadBC(const LoadBC& bc);
   LoadBC& operator=(const LoadBC& bc);

  }; // end class

} // end namespace
#endif

