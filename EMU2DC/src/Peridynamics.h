#ifndef EMU2DC_PERIDYNAMICS_H
#define EMU2DC_PERIDYNAMICS_H

#include <GlobalFlags.h>

namespace Emu2DC {

  typedef std::multimap<Node*,Node*> NodeFamily;
  typedef NodeFamily::const_iterator constNodeFamilyIterator;
  typedef NodeFamily::iterator NodeFamilyIterator;

  typedef std::multimap<Node*,Bond*> BondFamily;
  typedef BondFamily::const_iterator constBondFamilyIterator;
  typedef BondFamily::iterator BondFamilyIterator;
  typedef std::pair<Node*,Bond*> NodeBondPair;

  typedef std::vector<Bond*> BondArray;
  typedef std::vector<Bond*>::iterator BondIterator;

  class Peridynamics {

  public:
    Peridynamics(const GlobalFlags* flags);
    ~Peridynamics();

    void updateDisplacementVelocityVerlet();

  protected:

    void computeNodeFamily();
    void computeBondFamily();
    void getFamilyNodes(const Node* node,
		        NodeArray& familyNodes) const;
    void getFamilyBonds(const Node* node,
		        BondArray& familyBonds) const;
    void computeInternalForce(NodeArray& nodes);
    void computeBondForce(Bond* bond, 
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

    bool d_modified_mesh;
    int d_num_broken_bonds;
    int d_use_canonical_micromodulus;

    GlobalFlags* d_flags;

    // Keep material properties here for the time being
    Array3 d_damping;

    // Keep the node family here for the time being
    NodeFamily d_node_family;
    BondFamily d_bond_family;

  }; // end class
} // end namespace

#endif

