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

#ifndef __VAANGO_PERIDYNAMICS_BOND_DAMAGE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_BOND_DAMAGE__MODELH__

#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>

#include <vector>

namespace Vaango {

  class BondDamageModel: public PeridynamicsFailureModel {

  public:
         
    BondDamageModel(Uintah::ProblemSpecP& ps,
                    PeridynamicsFlags* flags);
    BondDamageModel(const BondDamageModel* cm);
    virtual ~BondDamageModel();

    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);
         
    /*! Initial computes and requires for the failure model */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    void initialize(const Uintah::Patch* patch,
                    const PeridynamicsMaterial* matl,
                    Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task  */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    // Basic failure model calculations
    void updateDamage(const Uintah::PatchSubset* patches,
                      const PeridynamicsMaterial* matl,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::DataWarehouse* new_dw);
                                     
    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

    BondDamageModel* clone();

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_BOND_DAMAGE_MODEL_H__

