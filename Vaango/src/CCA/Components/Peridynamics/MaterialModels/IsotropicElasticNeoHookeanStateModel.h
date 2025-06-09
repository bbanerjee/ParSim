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

#ifndef __VAANGO_PERIDYNAMICS_ISOTROPIC_ELASTIC_NEO_HOOKEAN_STATE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_ISOTROPIC_ELASTIC_NEO_HOOKEAN_STATE_MODEL_H__

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>

#include <Core/Grid/Variables/ComputeSet.h>    // For PatchSubset

/*!  The is essentially an isotropic compressible Neo-Hookean material model */

namespace Vaango {

  class IsotropicElasticNeoHookeanStateModel : public PeridynamicsMaterialModel {

  public:

    struct CMData {
      double bulkModulus;
      double shearModulus;
    };

  private:

    CMData d_cm;  // Constitutive model data

    // May be needed at a later stage (BB: 5/15/14)
    // friend const TypeDescription* fun_getTypeDescription(CMData* cm);

    IsotropicElasticNeoHookeanStateModel& operator=(const IsotropicElasticNeoHookeanStateModel& cm);

  public:
         
    // Default constructor
    IsotropicElasticNeoHookeanStateModel(Uintah::ProblemSpecP& ps,
                           PeridynamicsFlags* flags);

    // Copy constructor
    IsotropicElasticNeoHookeanStateModel(const IsotropicElasticNeoHookeanStateModel* cm);
    IsotropicElasticNeoHookeanStateModel(const IsotropicElasticNeoHookeanStateModel& cm) = default;

    // Make a clone of the constitutive model
    IsotropicElasticNeoHookeanStateModel* clone();

    virtual ~IsotropicElasticNeoHookeanStateModel();

    /*!  Output the problem spec for restart */
    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);

    /*! Initial computes and requires for the constitutive model */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    void initialize(const Uintah::Patch* patch,
                    const PeridynamicsMaterial* matl,
                    Uintah::DataWarehouse* new_dw);

    /*! Compute a stable initial timestep */
    void computeStableTimestep(const Uintah::Patch* patch,
                               const PeridynamicsMaterial* matl,
                               Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    /*! Compute the stress tensor */
    void computeStressTensor(const Uintah::PatchSubset* patches,
                             const PeridynamicsMaterial* matl,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw);

    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

  };
} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_ISOTROPIC_ELASTIC_NEO_HOOKEAN_STATE_MODEL_H__

