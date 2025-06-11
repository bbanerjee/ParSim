/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

// FrictionContactLRGuilkey.h

#ifndef __CAA_COMPONENTS_MPM_CONTACT_FRICTIONCONTACTLR_GUILKEY_H__
#define __CAA_COMPONENTS_MPM_CONTACT_FRICTIONCONTACTLR_GUILKEY_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>

#include <CCA/Ports/DataWarehouseP.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/**************************************

CLASS
   FrictionContactLRGuilkey

   This version of contact is based on John Nairn and Chad
   Hammerquist's 2019 manuscript that describes the use of logistic
   regression to find a common normal between objects, and uses
   particle geometry to find the most prominent portion of a particle
   at each node, and applies contact when contacting materials'
   prominences overlap.

GENERAL INFORMATION

   FrictionContactLRGuilkey.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   Contact_Model_Friction

DESCRIPTION
  One of the derived Contact classes.  This particular
  version is used to apply Coulombic frictional contact.

WARNING

****************************************/

class FrictionContactLRGuilkey : public Contact
{
  using Color = double;
  using Mu    = double;

private:
  // Coefficient of friction vs Color data
  std::vector<Color> d_color;
  std::vector<Mu> d_mu;
  std::vector<std::pair<Color, Mu>> d_color_mu;
  int d_material;

public:
  // Constructor
  FrictionContactLRGuilkey(const ProcessorGroup* myworld,
                           const MaterialManagerP& d_mat_manager,
                           const MPMLabel* labels,
                           const MPMFlags* flags,
                           ProblemSpecP& ps);

  // Destructor
  virtual ~FrictionContactLRGuilkey() = default;

  // Prevent copying/move of this class
  FrictionContactLRGuilkey(const FrictionContactLRGuilkey& con) = delete;
  FrictionContactLRGuilkey(FrictionContactLRGuilkey&& con)      = delete;
  FrictionContactLRGuilkey&
  operator=(const FrictionContactLRGuilkey& con) = delete;
  FrictionContactLRGuilkey&
  operator=(FrictionContactLRGuilkey&& con) = delete;

  virtual void
  setContactMaterialAttributes() override;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) override;

  // Basic contact methods
  virtual void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

  virtual void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;

private:
  virtual void
  exMomInterpolated(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  virtual void
  exMomIntegrated(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  inline double
  findMuFromColor(double color)
  {
    int n_entries = static_cast<int>(d_color.size());
    if (color >= d_color[n_entries - 1]) {
      return d_mu[n_entries - 1];
    }

    for (int ii = 1; ii < n_entries; ++ii) {
      if (color <= d_color[ii]) {
        double s = (d_color[ii] - color) / (d_color[ii] - d_color[ii - 1]);
        return d_mu[ii - 1] * s + d_mu[ii] * (1.0 - s);
      }
    }

    return d_mu[0];
  }
};
} // End namespace Uintah

#endif /* __CAA_COMPONENTS_MPM_CONTACT_FRICTIONCONTACTLR_GUILKEY_H__ */
