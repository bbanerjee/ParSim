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

// Friction.h

#ifndef __FRICTION_H__
#define __FRICTION_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

namespace Uintah {

class FrictionContact : public Contact
{
public:
  // Constructor
  FrictionContact(const ProcessorGroup* myworld,
                  const MaterialManagerP& mat_manager,
                  const MPMLabel* labels,
                  const MPMFlags* flags,
                  ProblemSpecP& ps);

  // Destructor
  virtual ~FrictionContact() = default;

  FrictionContact(const FrictionContact& con) = delete;
  FrictionContact(FrictionContact&& con)      = delete;
  FrictionContact&
  operator=(const FrictionContact& con) = delete;
  FrictionContact&
  operator=(FrictionContact&& con) = delete;

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

private:
  // Coefficient of friction
  double d_mu;

  // For hardcoded normals
  enum class NormalCoordSystem
  {
    NONE        = 0,
    CYLINDRICAL = 1,
    SPHERICAL   = 2,
    CARTESIAN   = 3
  };

  bool d_hardcodedNormals;
  std::vector<int> d_matIndex;

  // **TODO** Use map (key, value) instead
  std::vector<std::string> d_type;
  std::vector<NormalCoordSystem> d_coordType;
  std::vector<Point> d_center;
  std::vector<Vector> d_axisDir;
};
} // End namespace Uintah

#endif /* __FRICTION_H__ */
