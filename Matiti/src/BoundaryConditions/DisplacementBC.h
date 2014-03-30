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

