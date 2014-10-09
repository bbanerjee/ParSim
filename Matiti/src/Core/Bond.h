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

#ifndef MATITI_BOND_H
#define MATITI_BOND_H

#include <Pointers/NodeP.h>
#include <Pointers/WoodSP.h>
//#include <Containers/WoodSPArray.h>
#include <Woods/Wood.h>
#include <Containers/MaterialSPArray.h>
#include <Pointers/MaterialUP.h>
#include <Pointers/MaterialSP.h>
#include <MaterialModels/Material.h>
#include <Geometry/Vector3D.h>

#include <iostream>

namespace Matiti {

  class Bond 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Bond& bond);

  public: 

    Bond();
    Bond(const NodeP node1, const NodeP node2);
    Bond(const NodeP node1, const NodeP node2, const Material* mat);
    Bond(const NodeP node1, const NodeP node2, const Material* mat1, const Material* mat2);
    virtual ~Bond();

    /**
     * Compute volume weighted bond force
     */
    void computeInternalForce();
    void computeInternalForce(const MaterialSPArray& matList, const Vector3D& gridSize);

    /**
     * Compute volume weighted strain energy
     */
    double computeStrainEnergy() const;

    /**
     * Compute volume weighted micromodulus
     */
    double computeMicroModulus() const;

    /**
     * Compute the critical strain in the bond and flag if broken
     */
    bool checkAndFlagBrokenBond();
    bool checkAndFlagBrokenBond(const MaterialSPArray& matList); 

    /**
     * Set methods
     */
    void first(const NodeP node) { d_node1 = node; }
    void second(const NodeP node) { d_node2 = node; }
//    void material(MaterialUP& mat) { d_mat = std::move(mat); }
//    void material(const Material* mat) { d_mat->clone(mat); }
    void internalForce(const Vector3D& force)  { d_force = force; }
    void isBroken(bool broken) { d_broken = broken;}

    /**
     * Get methods
     */
    const NodeP first() const { return d_node1; }
//    const NodeP second() const { return d_node2; }
    const NodeP second() const { return d_node2; }
//    const Material* material() const {return d_mat.get(); }
    const MaterialSP materialS() const {return d_mat; }
    const Vector3D& internalForce() const { return d_force; }
    bool isBroken() const { return d_broken; }
    
    /**
     * Check if two bonds are identical
     *   True if the start and end points are the same 
     */
    bool operator==(const Bond& bond) const;

  private:
    NodeP d_node1;
    NodeP d_node2;
//    MaterialUP d_mat;
    MaterialSP d_mat;
//    WoodSP d_wood;
    Vector3D d_force;
    bool d_broken;

  };  // end class

} // end namespace

#endif
