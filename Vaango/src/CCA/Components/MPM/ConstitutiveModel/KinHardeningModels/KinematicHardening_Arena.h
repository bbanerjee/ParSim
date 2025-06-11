/*
 * The MIT License
 *
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

#ifndef __ARENA_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__
#define __ARENA_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardeningModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/AirEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/GraniteEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/WaterEOS.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class KinematicHardening_Arena
  \brief Backstress model for the pore pressure (Mason sand)
  \author Biswajit Banerjee,

  See ONR-MURI Reports from August 2015 through February 2016.
*/
/////////////////////////////////////////////////////////////////////////////

class KinematicHardening_Arena : public KinematicHardeningModel
{

  friend class InternalVar_Arena;

private:
  static const Uintah::Matrix3 Identity;

  struct CMData
  {
    double fluid_pressure_initial;
  };

  CMData d_cm;

  /* Tangent bulk modulus models for air, water, granite */
  AirEOS d_air;
  WaterEOS d_water;
  GraniteEOS d_granite;

  // Prevent copying of this class
  // KinematicHardening_Arena(const KinematicHardening_Arena &cm);
  KinematicHardening_Arena& operator=(const KinematicHardening_Arena& cm);

public:
  // constructors
  KinematicHardening_Arena(Uintah::ProblemSpecP& ps,
                           InternalVariableModel* intvar);
  KinematicHardening_Arena(const KinematicHardening_Arena* cm);

  // destructor
  ~KinematicHardening_Arena() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["Pf0"] = d_cm.fluid_pressure_initial;
    return params;
  }

  //////////
  /*! \brief Calculate the back stress */
  //////////
  void computeBackStress([[maybe_unused]] const ModelStateBase* state, [[maybe_unused]] const double& delT,
                         [[maybe_unused]] const Uintah::particleIndex idx,
                         [[maybe_unused]] const double& delLambda,
                         [[maybe_unused]] const Uintah::Matrix3& df_dsigma_new,
                         [[maybe_unused]] const Uintah::Matrix3& backStress_old,
                         [[maybe_unused]] Uintah::Matrix3& backStress_new) override
  {
  }

  void computeBackStress(const ModelStateBase* state,
                         Uintah::Matrix3& backStress_new) override;

  void eval_h_beta([[maybe_unused]] const Uintah::Matrix3& df_dsigma,
                   [[maybe_unused]] const ModelStateBase* state,
                   [[maybe_unused]] Uintah::Matrix3& h_beta) override
  {
  }

private:
  //------------------------------------------------------
  // Newton solve for pressure
  //------------------------------------------------------
  double computePressureUnloaded(const double& ev_p, const double& S0,
                                 const double& phi0, const double& p0,
                                 const double& p_init);
  double computeGpByDgp(const double& p, const double& ev_p, const double& S0,
                        const double& phi0, const double& p0);

  //------------------------------------------------------
  // Bisection solve for pressure
  //------------------------------------------------------
  double computePressureBisection(const double& ev_p, const double& S0,
                                  const double& phi0, const double& p0,
                                  const double& p_init);

  double computeGp(const double& p, const double& ev_p, const double& S0,
                   const double& phi0, const double& p0);

public:
  // We use the Matrix3 pBackStressLabel instead of pZetaLabel.
  // pBackStressLabel is defined in the base class.
  // const Uintah::VarLabel*   pZetaLabel;
  // const Uintah::VarLabel*   pZetaLabel_preReloc;

  // Add particle state for these labels
  void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                        std::vector<const Uintah::VarLabel*>& to) override
  {
    from.push_back(pBackStressLabel);
    to.push_back(pBackStressLabel_preReloc);
  }

  /**
   * Initialize pore pressure label
   */
  void initializeLocalMPMLabels()
  {
    pBackStressLabel = Uintah::VarLabel::create(
      "p.porePressure",
      Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
    pBackStressLabel_preReloc = Uintah::VarLabel::create(
      "p.porePressure+",
      Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  }

  /**
   * Set up task graph for initialization
   */
  void addInitialComputesAndRequires(
    Uintah::Task* task, const Uintah::MPMMaterial* matl,
    [[maybe_unused]] const Uintah::PatchSet* patch) const override
  {
    const Uintah::MaterialSubset* matlset = matl->thisMaterial();
    task->computes(pBackStressLabel, matlset);
  }

  /**
   *  Actually initialize the pore pressure
   */
  void initializeLocalVariables(
    [[maybe_unused]] const Uintah::Patch* patch, Uintah::ParticleSubset* pset,
    Uintah::DataWarehouse* new_dw,
    [[maybe_unused]] Uintah::constParticleVariable<double>& pVolume) override
  {
    Uintah::ParticleVariable<Uintah::Matrix3> pBackStress;
    new_dw->allocateAndPut(pBackStress, pBackStressLabel, pset);

    for (int idx : *pset) {
      pBackStress[idx] = -d_cm.fluid_pressure_initial * Identity;
    }
  }

  /**
   * Set up task graph for parameter copying to new datawarehouse
   */
  void addComputesAndRequires(Uintah::Task* task,
                              const Uintah::MPMMaterial* matl,
                              [[maybe_unused]] const Uintah::PatchSet* patches) const override
  {
    const Uintah::MaterialSubset* matlset = matl->thisMaterial();
    task->needs(Uintah::Task::OldDW, pBackStressLabel, matlset,
                   Uintah::Ghost::None);
    task->computes(pBackStressLabel_preReloc, matlset);
  }

  /**
   *  Update the pore pressure
   */
  void computeBackStress([[maybe_unused]] Uintah::ParticleSubset* pset,
                         [[maybe_unused]] Uintah::DataWarehouse* old_dw,
                         [[maybe_unused]] Uintah::DataWarehouse* new_dw)
  {
  }

  void getBackStress(
    Uintah::ParticleSubset* pset, Uintah::DataWarehouse* old_dw,
    Uintah::constParticleVariable<Uintah::Matrix3>& pBackStress) override
  {
    old_dw->get(pBackStress, pBackStressLabel, pset);
  }

  void allocateAndPutBackStress(
    Uintah::ParticleSubset* pset, Uintah::DataWarehouse* new_dw,
    Uintah::ParticleVariable<Uintah::Matrix3>& pBackStress_new) override
  {
    new_dw->allocateAndPut(pBackStress_new, pBackStressLabel_preReloc, pset);
  }
};

} // End namespace Uintah

#endif // __ARENA_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__
