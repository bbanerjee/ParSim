#ifndef EMU2DC_BODY_H
#define EMU2DC_BODY_H

#include <Domain.h>
#include <SimulationState.h>
#include <FamilyComputer.h>
#include <Material.h>
#include <MaterialSPArray.h>
#include <CrackSPArray.h>
#include <NodePArray.h>
#include <ElementPArray.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <iostream>
#include <map>

namespace Emu2DC {

  class Body
  {
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Body& body);

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

    void removeBondsIntersectedByCracks();

    inline int id() const {return d_id;}
    inline void id(const int& id) {d_id = id;}

    // **WARNING** One mat for now.  A body can have more than one material. Also the
    // materials can be PerMaterial, MPMMaterial, or RigidMaterial.
    inline int matID() const {return d_mat_id;}
    const NodePArray& nodes() const {return d_nodes;}
    const ElementPArray& elements() const {return d_elements;}
    const FamilyComputer& familyComputer() const {return d_family_computer;}
    const Vector3D& initialVelocity() const {return d_initial_velocity;}
    const Vector3D& bodyForce() const {return d_body_force;}
    const CrackSPArray& cracks() const {return d_cracks;}
   

  protected:

    void readNodeFile(const std::string& fileName, const int dim);
    void setInitialNodeHorizon(const double horizon);
    void assignNodeMaterial(const MaterialSPArray& matList);

    void readElementFile(const std::string& fileName);
    void computeNodalVolumes();
    void initializeFamilyComputer(const Domain& domain);

  private:

    int d_id;
    int d_mat_id;
    NodePArray d_nodes;
    ElementPArray d_elements;

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

    FamilyComputer d_family_computer;

    Vector3D d_initial_velocity; // Initial velocity
    Vector3D d_body_force;       // Gravity (essentially)

    CrackSPArray d_cracks;

  };
} // end namespace

#endif
