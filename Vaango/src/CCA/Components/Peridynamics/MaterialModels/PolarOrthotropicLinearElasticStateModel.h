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

#ifndef __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>

#include <Core/Grid/Variables/ComputeSet.h>    // For PatchSubset
#include <Core/Math/SymmMatrix6.h>

/*!  A polar orthotropic linear elastic material model which uses a Green Naghdi stress rate */

namespace Vaango {

  class PolarOrthotropicLinearElasticStateModel : public PeridynamicsMaterialModel {

  public:

    struct CMData {
      Uintah::Point top;
      Uintah::Point bottom;
      double Er;
      double Etheta;
      double Ez;
      double nuthetar;
      double nuzr;
      double nuztheta;
      double Gthetaz;
      double Gzr;
      double Grtheta;
      Uintah::SymmMatrix6 stiffnessMatrix;
    };

  private:

    CMData d_cm;  // Constitutive model data

    // Prevent assignment
    PolarOrthotropicLinearElasticStateModel& operator=(const PolarOrthotropicLinearElasticStateModel& cm);

  public:
         
    // Default constructor
    PolarOrthotropicLinearElasticStateModel(Uintah::ProblemSpecP& ps,
                                            PeridynamicsFlags* flags);

    // Copy constructor
    PolarOrthotropicLinearElasticStateModel(const PolarOrthotropicLinearElasticStateModel* cm);

    // Make a clone of the constitutive model
    PolarOrthotropicLinearElasticStateModel* clone();

    virtual ~PolarOrthotropicLinearElasticStateModel();

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
      

#endif  // __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__

