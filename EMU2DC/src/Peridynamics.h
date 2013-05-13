#ifndef EMU2DC_PERIDYNAMICS_H
#define EMU2DC_PERIDYNAMICS_H

#include <SimulationState.h>
#include <Time.h>
#include <Output.h>
#include <Domain.h>
#include <MaterialSPArray.h>
#include <BodySPArray.h>
#include <HorizonComputer.h>

#include <PeridynamicsTypes.h>
#include <NodeP.h>
#include <NodePArray.h>
#include <BondP.h>
#include <BondPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Emu2DC {

  class Peridynamics {

  public:
    Peridynamics();
    ~Peridynamics();

    void problemSetup(Uintah::ProblemSpecP& ps);

    void run();

    void updateDisplacementVelocityVerlet();

  protected:

    void computeInternalForce(const NodePArray& nodeList);

    void computeNodeFamily();
    void computeBondFamily();
    void getFamilyNodes(const NodeP node,
		        NodePArray& familyNodes) const;
    void getFamilyBonds(const NodeP node,
		        BondPArray& familyBonds) const;
    void computeBondForce(BondP bond, 
		          Array3& bondForce,
			  double& bondLengthInit,
			  double& bondLengthNew,
			  double& bondStrain,
			  double& bondStrainEnergy,
			  double& micromodulus);
    double computeMicromodulus(const double& bondLengthInitial, 
		               const double& horizonRadius,
			       const double& youngsModulus);

    void integrateNodalAcceleration();
    void breakBonds();

  private:

    Time d_time;
    Output d_output;
    SimulationState d_state;
    Domain d_domain;
    MaterialSPArray d_mat_list;
    BodySPArray d_body_list;

    int d_num_broken_bonds;

    // Keep the node family here for the time being
    NodeFamily d_node_family;
    BondFamily d_bond_family;

  }; // end class
} // end namespace

#endif

