#ifndef MATITI_BODY_H
#define MATITI_BODY_H

#include <Domain.h>
#include <SimulationState.h>
//#include <FamilyComputer.h>
#include <Material.h>
#include <MaterialSPArray.h>
#include <NodePArray.h>
//#include <ElementPArray.h>
#include <InitialConditions.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <iostream>
#include <map>

namespace BrMPM {

  class Body
  {
  public:  

    friend std::ostream& operator<<(std::ostream& out, const BrMPM::Body& body);

  public:
   
    Body();
    virtual ~Body();

    void initialize(Uintah::ProblemSpecP& ps,
                    const Domain& domain,
                    const SimulationState& state, 
                    const MaterialSPArray& matList);

    void createInitialFamily(const Domain& domain);
    void updateFamily(const Domain& domain);
    void printFamily();

    inline int id() const {return d_id;}
    inline void id(const int& id) {d_id = id;}

    // **WARNING** One mat for now.  A body can have more than one material. Also the
    // materials can be PerMaterial, MPMMaterial, or RigidMaterial.
    inline int matID() const {return d_mat_id;}
    const NodePArray& nodes() const {return d_nodes;}
    const ElementPArray& elements() const {return d_elements;}
    const FamilyComputer& familyComputer() const {return d_family_computer;}

    /**
     * Get methods for initial conditions
     */
    const Vector3D& initialVelocity() const {return d_ic.initialVelocity();}
    const Vector3D& bodyForce() const {return d_ic.bodyForce();}

    /**
     * Compute methods for initial bond removal
     */
    void removeBondsIntersectedByCracks(){d_ic.removeBondsIntersectedByCracks(d_nodes);}

     /**
      * Update damage index of all nodes in the body
      */
    void updateDamageIndex() const;

  protected:

    void readMaterialInput(Uintah::ProblemSpecP& ps,
                           const MaterialSPArray& matList);

    void readNodeFile(const std::string& fileName);
    void setInitialNodeHorizon(const double horizon);
    void assignNodeMaterial(const MaterialSPArray& matList);

    void readElementFile(const std::string& fileName);
    void computeNodalVolumes();
    void initializeFamilyComputer(const Domain& domain);

  private:

    int d_id;
    int d_mat_id;
    std::string d_mat_dist;
    double d_mat_cov;
    double d_mat_seed;

    NodePArray d_nodes;
    ElementPArray d_elements;

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

    FamilyComputer d_family_computer;

    InitialConditions d_ic;

  };
} // end namespace

#endif
