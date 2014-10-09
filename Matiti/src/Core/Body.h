/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MATITI_BODY_H
#define MATITI_BODY_H

#include <Core/Domain.h>
#include <Core/SimulationState.h>
#include <Core/FamilyComputer.h>
#include <MaterialModels/Material.h>
#include <MaterialModels/Density.h>
#include <Containers/MaterialSPArray.h>
#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>
#include <Containers/LoadBCSPArray.h>
#include <Containers/DispBCSPArray.h>
#include <Containers/WoodSPArray.h>
#include <Pointers/DensitySP.h>
#include <Pointers/WoodSP.h>
#include <BoundaryConditions/InitialConditions.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <fstream>
#include <iostream>
#include <map>

namespace Matiti {

  class Body
  {
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Body& body);

  public:
   
    Body();
    virtual ~Body();

    void initialize(Uintah::ProblemSpecP& ps,
                    const Domain& domain,
                    const SimulationState& state, 
                    const MaterialSPArray& matList);

    void initialize(int materialId,
                    NodePArray& nodes,
                    ElementPArray& elements,
                    const Vector3D& gridSize,
                    const Domain& domain,
                    const SimulationState& state, 
                    const MaterialSPArray& matList,
                    const InitialConditions& ic);

    void createInitialFamily(const Domain& domain);
    void updateFamily(const Domain& domain);
    void printFamily();

    void applyDisplacementBC();

    inline int id() const {return d_id;}
    inline void id(const int& id) {d_id = id;}

    // **WARNING** One mat for now.  A body can have more than one material. Also the
    // materials can be PerMaterial, MPMMaterial, or RigidMaterial.
    inline int matID() const {return d_mat_id;}

    const Vector3D& gridSize() const {return d_grid_size;}
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

    void setInitialNodeHorizon(const double horizon);
    void assignNodeMaterial(const MaterialSPArray& matList);

    void computeNodalVolumes();
    void computeNodalDensity(const MaterialSPArray& matList);
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
    LoadBCSPArray d_load_bcs;
    DispBCSPArray d_disp_bcs;

    Vector3D d_grid_size;

  };
} // end namespace

#endif
