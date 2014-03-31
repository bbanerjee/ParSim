#ifndef MATITI_PERIDYNAMICS_H
#define MATITI_PERIDYNAMICS_H

#include <Core/SimulationState.h>
#include <Core/Time.h>
#include <InputOutput/OutputVTK.h>
#include <Core/Domain.h>
#include <Containers/MaterialSPArray.h>
#include <Containers/BodySPArray.h>
#include <Core/HorizonComputer.h>
#include <BoundaryConditions/VelocityBC.h>

#include <Types/PeridynamicsTypes.h>
#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Pointers/BondP.h>
#include <Containers/BondPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {

  class Peridynamics {

  public:
    Peridynamics();
    ~Peridynamics();

    void problemSetup(Uintah::ProblemSpecP& ps);

    void run();

  protected:

    void applyInitialConditions();

    void computeInternalForce(const NodeP& node,
                              Vector3D& internalForce);

    void integrateNodalAcceleration(const Vector3D& velocityOld,
                                    const Vector3D& accelerationOld,
                                    double delT,
                                    Vector3D& velocityNew);

    void integrateNodalVelocity(const Vector3D& displacementOld,
                                const Vector3D& velocityOld,
                                double delT,
                                Vector3D& displacementNew);

    double computeMicromodulus(const double& bondLengthInitial, 
		               const double& horizonRadius,
			       const double& youngsModulus);

    void breakBonds(const NodePArray& nodes);

    void checkMemoryUsage(double& resident_mem, double& shared_mem);

  private:

    Time d_time;
    OutputVTK d_output;
    SimulationState d_state;
    Domain d_domain;
    MaterialSPArray d_mat_list;
    BodySPArray d_body_list;
    VelocityBC d_velocitybc;

    int d_num_broken_bonds;

  }; // end class
} // end namespace

#endif

