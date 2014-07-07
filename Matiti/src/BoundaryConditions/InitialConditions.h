#ifndef MATITI_INITIALCONDITIONS_H
#define MATITI_INITIALCONDITIONS_H

#include <Containers/NodePArray.h>
#include <Containers/CrackSPArray.h>
#include <Geometry/Vector3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {
  
  class InitialConditions {
  
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::InitialConditions& ic);

  public:

    InitialConditions();
    InitialConditions(const Vector3D& initialVel,
                      const Vector3D& gravity);
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

    /**
     * Set methods
     */
    void initialVelocity(const Vector3D& vel) {d_initial_velocity = vel;}
    void bodyForce(const Vector3D& gravity) {d_body_force = gravity;}

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

