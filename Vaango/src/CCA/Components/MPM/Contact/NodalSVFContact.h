/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

// NodalSVFContact.h

#ifndef __Nodal_SVF_H__
#define __Nodal_SVF_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {
/**************************************

CLASS
   NodalSVFContact

 // This is a contact model developed by Peter Mackenzie. Details about the
 // derivation can be found in "Modeling Strategies for Multiphase Drag
 // Interactions Using the Material Point Method" (Mackenzie, et al; 2011).
 // Of the interaction models proposed in their paper, this particular
 // contact model for MPM in Uintah can simulate the Nodal Bang Bang method
 // OR the Nodal Smoothed Volume (SVF) Fraction Method. The Nodal Bang Bang
 // method is less expensive than Nodal SVF, but much less accurate. These
 // two methods, which are of the Node-based type, register interaction
 // proportional to the cell volume (dx*dy*dz) between two or more phases.
 // As a result, over-estimates of the interaction occur in cells where the
 // interacting materials do not completely occupy the computational cell.
 // Interaction in this particular code module is quantified by an interaction
 // parameter, mu (N/m^4), and the velocity difference between the phases. Other
 // simple Coulomb friction or viscous fluid interaction models can be
 // substituted.

GENERAL INFORMATION

   NodalSVFContact.h

   Coded by Cherie L. Kunkel
   Department of Mechanical Engineering
   University of Utah

   Based on the Model Developed by Peter Mackenzie

KEYWORDS
   Contact_Model_Nodal_SVF_Smoothed_Volume_Fraction

WARNING

****************************************/

class NodalSVFContact : public Contact
{
private:

  // PARAMETERS UNIQUE TO THIS MODEL FROM UPS FILE
  double d_myu{0.0};
  bool d_svf{false};
  std::vector<int> d_materials;

public:
  // Constructor
  NodalSVFContact(const ProcessorGroup* myworld,
                  const MaterialManagerP& mat_manager,
                  const MPMLabel* labels,
                  const MPMFlags* flags,
                  ProblemSpecP& ps);

  // Destructor
  virtual ~NodalSVFContact() = default;

  // Prevent copying/move of this class
  NodalSVFContact(const NodalSVFContact& con) = delete;
  NodalSVFContact(NodalSVFContact&& con)      = delete;
  NodalSVFContact&
  operator=(const NodalSVFContact& con) = delete;
  NodalSVFContact&
  operator=(NodalSVFContact&& con) = delete;

  virtual void
  setContactMaterialAttributes() override;

  void
  outputProblemSpec(ProblemSpecP& ps) override;

  void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

  void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;
};
} // End namespace Uintah

#endif /* __Nodal_SVF_H__ */
