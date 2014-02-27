#ifndef __MATITI_TRACTION_BC_H__
#define __MATITI_TRACTION_BC_H__

#include <BoundaryConditions/LoadBC.h>

namespace Matiti {
  
  class TractionBC : public LoadBC {
  
  public:

    TractionBC();
    ~TractionBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elems);

  protected:

    void computeExtForceDensity(const SCIRun::Vector& extForce,
                                NodePArray& surfaceNodes, 
                                ElementPArray& elems);

  private:

   // prevent copying
   TractionBC(const TractionBC& dyna);
   TractionBC& operator=(const TractionBC& dyna);

  }; // end class

} // end namespace
#endif

