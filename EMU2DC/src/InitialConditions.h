#ifndef EMU2DC_INITIALCONDITIONS_H
#define EMU2DC_INITIALCONDITIONS_H

#include <NodePArray.h>
#include <CrackSPArray.h>
#include <Geometry/Vector3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Emu2DC {
  
  class InitialConditions {
  
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::InitialConditions& ic);

  public:

    InitialConditions();
    ~InitialConditions();

    void readInitialConditions(Uintah::ProblemSpecP& ps);

    void applyInitialVelocity(NodePArray& nodes);
    void applyBodyForce(NodePArray& nodes);
    void removeBondsIntersectedByCracks(NodePArray& nodes);

    /**
     * Get methods
     */
    const Vector3D& initialVelocity() const {return d_initial_velocity;}
    const Vector3D& bodyForce() const {return d_body_force;}
    const CrackSPArray& cracks() const {return d_cracks;}

  private:

    Vector3D d_initial_velocity; // Initial velocity
    Vector3D d_body_force;       // Gravity (essentially)
    CrackSPArray d_cracks;       // Initial crack geometries

    // prevent copying
    InitialConditions(const InitialConditions& bc);
    InitialConditions& operator=(const InitialConditions& bc);

  }; // end class

} // end namespace
#endif

