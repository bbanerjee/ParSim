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

#ifndef __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class Task;
  class Patch;
  class VarLabel;
  class DataWarehouse;
  class ParticleSubset;
  class ParticleVariableBase;
}

namespace Vaango {

  class PeridynamicsLabel;
  class PeridynamicsFlags;
  class PeridynamicsMaterial;

  class PeridynamicsDamageModel {

  public:
         
    /*! Default constructor */
    PeridynamicsDamageModel(PeridynamicsLabel* labels, PeridynamicsFlags* flags);

    /*! Copy constructor */
    PeridynamicsDamageModel(const PeridynamicsDamageModel* cm);

    /*! Make a clone of the constitutive model */
    virtual PeridynamicsDamageModel* clone() = 0;

    /*! Destructor */
    virtual ~PeridynamicsDamageModel();

    /*! Write out the input data for restart purposes */
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps,
                                   bool output_cm_tag = true) = 0;

    /*! Set up initial computes and requires for the damage model compute task */
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const PeridynamicsMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    /*! Actually initialize the variables used in the damage models */
    virtual void initialize(const Uintah::Patch* patch,
                            const PeridynamicsMaterial* matl,
                            Uintah::DataWarehouse* new_dw) = 0;

    /*! Set up the computes and requires for the task that computes the
      damage tensor */
    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const PeridynamicsMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    /*! Actually compute damage tensor */
    virtual void computeDamageTensor(const Uintah::PatchSubset* patches,
                                     const PeridynamicsMaterial* matl,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw) = 0;
                                     
    /*! Tag variables that are to be relocated */
    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to) = 0;


  protected:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flags;

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__

