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

#include <CCA/Components/Peridynamics/Core/PeridynamicsDomainBoundCond.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <vector>
#include <iostream>

using namespace Vaango;

using Uintah::Patch;
using Uintah::NCVariable;
using Uintah::Iterator;
using Uintah::NodeIterator;
using Uintah::BoundCondBase;
using Uintah::BoundCondBaseSP;
using Uintah::BoundCond;

using Uintah::Vector;
using Uintah::IntVector;

PeridynamicsDomainBoundCond::PeridynamicsDomainBoundCond()
{
}

PeridynamicsDomainBoundCond::~PeridynamicsDomainBoundCond()
{
}

/*-------------------------------------------------------------------------------------------
 * setBoundaryCondition: <Vector>
 *  Sets the boundary condition for node-centered Uintah::Vector variables based on bc_type 
 *  Two types of boundary conditions are allowed:
 *    Dirichlet:  "Velocity"  bcs
 *    Symmetry:   "Symmetry" bcs
 *  Note: Only linear interpolation allowed at this stage
 *-------------------------------------------------------------------------------------------
 */
void 
PeridynamicsDomainBoundCond::setBoundaryCondition(const Patch* patch,
                                                  int matlIndex,
                                                  const std::string& bc_type, 
                                                  NCVariable<Vector>& variable,
                                                  std::string interp_type)
{
  for (auto face = Patch::startFace; face <= Patch::endFace; 
            face = Patch::nextFace(face)) {

    IntVector oneCell = patch->faceDirection(face);

    if (patch->getBCType(face) == Patch::None) {

      int numChildren = patch->getBCDataArray(face)->getNumberChildren(matlIndex);

      for (int child = 0; child < numChildren; child++) {

        Iterator nbound_ptr, dummy;

        BoundCondBaseSP bcb = 
          patch->getArrayBCValues(face, matlIndex, bc_type, dummy, nbound_ptr, child);
        if (!bcb) continue;

        // Velocity BC
        if (bc_type == "Velocity") {

          // Cast the boundary condition into a vector type BC
          BoundCond<Vector>::BoundCondP bc = std::dynamic_pointer_cast<BoundCond<Vector> >(bcb); 
          if (bc) {
            if (bc->getBCType() == "Dirichlet") {
              for (nbound_ptr.reset();!nbound_ptr.done();nbound_ptr++){ 
                variable[*nbound_ptr] = bc->getValue();
              }
            } // end if Dirichlet BC

          } // end if (bc)

        } else if (bc_type == "Symmetric") {

          if (bcb->getBCType() == "symmetry") {

            if (face == Patch::xplus || face == Patch::xminus){
              for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                IntVector nd = *nbound_ptr;
                variable[nd] = Vector(0.0, variable[nd].y(), variable[nd].z());
              }
            } // xplus/xminus faces

            if (face == Patch::yplus || face == Patch::yminus){
              for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                IntVector nd = *nbound_ptr;
                variable[nd] = Vector(variable[nd].x(),0.,variable[nd].z());
              }
            }  // yplus/yminus faces

            if (face == Patch::zplus || face == Patch::zminus){
              for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++){
                IntVector nd = *nbound_ptr;
                variable[nd] = Vector(variable[nd].x(), variable[nd].y(),0.);
              }
            } // zplus/zminus
          } // end if bc_type_new = symmetry
        } // end if Symmetric BC
 
      } // end child loop
    } // end if Patch::None 
  } // end face loop
}

/*-------------------------------------------------------------------------------------------
 * setBoundaryCondition: <double>
 *  Sets the boundary condition for node-centered double variables based on bc_type 
 *  Two types of boundary conditions are allowed:
 *    Dirichlet:  "Velocity"  bcs
 *    Symmetry:   "Symmetry" bcs
 *  Note: Only linear interpolation allowed at this stage
 *-------------------------------------------------------------------------------------------
 */
void 
PeridynamicsDomainBoundCond::setBoundaryCondition(const Patch* patch, 
                                                  int matlIndex,
                                                  const std::string& bc_type, 
                                                  NCVariable<double>& variable,
                                                  std::string interp_type)
{
  // Loop through faces of the domain box
  for (auto face = Patch::startFace; face <= Patch::endFace; face = Patch::nextFace(face)) {

    IntVector oneCell = patch->faceDirection(face);

    if (patch->getBCType(face) == Patch::None) {

      int numChildren = patch->getBCDataArray(face)->getNumberChildren(matlIndex);

      for (int child = 0; child < numChildren; child++) {

        // Get the BC data array
        Iterator nbound_ptr, dummy;
        BoundCondBaseSP bcb = 
          patch->getArrayBCValues(face, matlIndex, bc_type, dummy, nbound_ptr, child);
        if (!bcb) continue;

        // Cast the boundary condition into a double type
        BoundCond<double>::BoundCondP bc = std::dynamic_pointer_cast<BoundCond<double> >(bcb);
        if (!bc) continue;

        if (bc_type=="Pressure" || bc_type=="Temperature") {

          if (bc->getBCType() == "Dirichlet") {
            double bcv = bc->getValue();
            for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
              variable[*nbound_ptr] = bcv;
            }
          } // end if Dirichlet
            
          if (bc->getBCType() == "Neumann"){
            Vector deltax = patch->dCell();
            double dx = -9;
            IntVector off(-9,-9,-9);
            if (face == Patch::xplus){
              dx = deltax.x();
              off=IntVector(1,0,0);
            } else if (face == Patch::xminus) {
              dx = deltax.x();
              off=IntVector(-1,0,0);
            } else if (face == Patch::yplus) {
              dx = deltax.y();
              off=IntVector(0,1,0);
            } else if (face == Patch::yminus) {
              dx = deltax.y();
              off=IntVector(0,-1,0);
            } else if (face == Patch::zplus) {
              dx = deltax.z();
              off=IntVector(0,0,1);
            } else if (face == Patch::zminus) {
              dx = deltax.z();
              off=IntVector(0,0,-1);
            }
              
            double gradv = bc->getValue();
            for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
	      IntVector nd = *nbound_ptr;
              variable[nd] = variable[nd-off] - gradv*dx;
	    }
          } // end if Neumann
            
        } // end if Pressure or Temperature bc

      }  // end child loop
    } // end if Patch::None 
  } // end face loop
}
