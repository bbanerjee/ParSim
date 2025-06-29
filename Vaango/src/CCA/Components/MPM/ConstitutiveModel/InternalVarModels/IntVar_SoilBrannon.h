/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __SOIL_MODEL_BRANNON_KAPPA_INT_VAR_MODEL_H__
#define __SOIL_MODEL_BRANNON_KAPPA_INT_VAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_SoilModelBrannon.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class IntVar_SoilBrannon
  \brief The evolution of the kappa hardening internal variable in the
         Arenisca model

  Reference: Arenisca manual.

  The rate of change of kappa is given by

   dkappa/deps_v = F(kappa) - G(eps_v) - H(eps_v)

  where

   eps_v = volumetric plastic strain

  and

   F(kappa) = 1.0/(p1*p3)*exp(-p1*kappa - p0)
   G(eps_v) = B1 exp(p3 + p4 + eps_v)/[exp(p3 + p4 + eps_v) - 1]^2
   H(eps_v) = B1 exp(p3 + eps_v)/[exp(p3 + eps_v) - 1]^2

  with
   B1 = 3 B0 [exp(p3 + p4) - 1]


  The incremental update of the consolidation pressure is given by

     kappa_{n+1} = kappa_n + [F(kappa_{n+1}) - G(eps_v) - H(eps_v)] Delta eps_v
*/
////////////////////////////////////////////////////////////////////////////

class IntVar_SoilBrannon : public InternalVariableModel
{

public:
  // Internal variables
  const Uintah::VarLabel* pKappaLabel;
  const Uintah::VarLabel* pKappaLabel_preReloc;

  // Return the internal variable labels
  std::vector<const Uintah::VarLabel*>
  getLabels() const override
  {
    std::vector<const Uintah::VarLabel*> labels;
    labels.push_back(pKappaLabel); // Branch point
    labels.push_back(pKappaLabel_preReloc);

    return labels;
  }

  // constructors
  IntVar_SoilBrannon(Uintah::ProblemSpecP& ps);
  IntVar_SoilBrannon(const IntVar_SoilBrannon* cm);
  IntVar_SoilBrannon&
  operator=(const IntVar_SoilBrannon& cm) = delete;

  // destructor
  ~IntVar_SoilBrannon() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["p0"]     = d_p0;
    params["p1"]     = d_p1;
    params["p3"]     = d_p3;
    params["p4"]     = d_p4;
    params["B0"]     = d_B0;
    params["CR"]     = d_Cr;
    params["FSLOPE"] = d_fSlope;
    params["PEAKI1"] = d_peakI1;
    return params;
  }

  void
  addParticleState(std::vector<const Uintah::VarLabel*>& from,
                   std::vector<const Uintah::VarLabel*>& to) override;

  // Computes and requires for internal evolution variables
  void
  addInitialComputesAndRequires(Uintah::Task* task,
                                const Uintah::MPMMaterial* matl,
                                const Uintah::PatchSet* patches) override;

  void
  initializeInternalVariable(Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw) override;
  void
  initializeInternalVariable([[maybe_unused]] const Uintah::Patch* patch,
                             [[maybe_unused]] const Uintah::MPMMaterial* matl,
                             [[maybe_unused]] Uintah::ParticleSubset* pset,
                             [[maybe_unused]] Uintah::DataWarehouse* new_dw,
                             [[maybe_unused]] Uintah::MPMLabel* lb,
                             [[maybe_unused]] ParamMap& params) override
  {
  }

  void
  addComputesAndRequires(Uintah::Task* task,
                         const Uintah::MPMMaterial* matl,
                         const Uintah::PatchSet* patches) override;

  /* Get one (possibly composite) internal variable */
  template <typename T>
  void
  getInternalVariable(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::constParticleVariable<T>& intvar);

  /* Get multiple local <int/double/Vector/Matrix3> internal variables */
  template <typename T>
  std::vector<Uintah::constParticleVariable<T>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw);

  /* Allocate one (possibly composite) internal variable */
  template <typename T>
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 Uintah::ParticleVariable<T>& intvar);

  /* Allocate multiple local <int/double/Vector/Matrix3> internal variables */
  template <typename T>
  void
  allocateAndPutInternalVariable(
    Uintah::ParticleSubset* pset,
    Uintah::DataWarehouse* new_dw,
    std::vector<Uintah::ParticleVariable<T>>& pVars);

  ///////////////////////////////////////////////////////////////////////////
  /*! \brief Compute the internal variable */
  template <typename T>
  void
  evolveInternalVariable(Uintah::particleIndex pidx,
                         const ModelStateBase* state,
                         Uintah::constParticleVariable<T>& var_old,
                         Uintah::ParticleVariable<T>& var);
  double
  computeInternalVariable(const std::string& label,
                          const ModelStateBase* state_old,
                          const ModelStateBase* state_cur) const override;

  ///////////////////////////////////////////////////////////////////////////
  // Compute derivative of internal variable with respect to volumetric
  // elastic strain
  double
  computeVolStrainDerivOfInternalVariable([[maybe_unused]] const std::string& label,
                                          [[maybe_unused]] const ModelStateBase*) const override
  {
    return 0.0;
  }

  void
  allocateCMDataAddRequires(Uintah::Task* task,
                            const Uintah::MPMMaterial* matl,
                            const Uintah::PatchSet* patch,
                            Uintah::MPMLabel* lb) override;

  void
  allocateCMDataAdd(
    Uintah::DataWarehouse* new_dw,
    Uintah::ParticleSubset* addset,
    std::map<const Uintah::VarLabel*, Uintah::ParticleVariableBase*>* newState,
    Uintah::ParticleSubset* delset,
    Uintah::DataWarehouse* old_dw) override;

  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleVariableBase& intvar) override;

  void
  allocateAndPutRigid([[maybe_unused]] Uintah::ParticleSubset* pset,
                      [[maybe_unused]] Uintah::DataWarehouse* new_dw,
                      [[maybe_unused]] Uintah::constParticleLabelVariableMap& intvars) override
  {
  }

private:
  // Model parameters
  double d_p0;
  double d_p1;
  double d_p3;
  double d_p4;
  double d_B0;
  double d_Cr;
  double d_fSlope;
  double d_peakI1;

  //--------------------------------------------------------------------------------------
  // Compute kappa_new from the function X1(kappa_{n+1})
  //  where
  //       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //
  // ** NOTE** (should be replaced with function pointers)
  //--------------------------------------------------------------------------------------
  double
  computeKappaFromX1(const double& kappa_old,
                     const double& epsv,
                     const double& deltaEpsv,
                     const double& tolerance,
                     const int& maxiter) const;

  //--------------------------------------------------------------------------------------
  // Compute kappa_new from the function X2(kappa_{n+1})
  //  where
  //       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //
  // ** NOTE** (should be replaced with function pointers)
  //--------------------------------------------------------------------------------------
  double
  computeKappaFromX2(const double& kappa_old,
                     const double& epsv,
                     const double& deltaEpsv,
                     const double& tolerance,
                     const int& maxiter) const;

  //--------------------------------------------------------------------------------------
  // Compute the function X1(kappa_{n+1})
  //  where
  //       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //--------------------------------------------------------------------------------------
  double
  computeX1(const double& kappa_old,
            const double& kappa_new,
            const double& G,
            const double& H,
            const double& delEpsv) const;

  //--------------------------------------------------------------------------------------
  // Compute the function dX1/dkappa(kappa_{n+1})
  //  where
  //       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //--------------------------------------------------------------------------------------
  double
  computeDerivX1dkappa(const double& kappa_old,
                       const double& kappa_new,
                       const double& delEpsv) const;

  //--------------------------------------------------------------------------------------
  // Compute the value of kappa at which function X1 is a minimum
  //  where
  //       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //--------------------------------------------------------------------------------------
  double
  computeKappaAtX1Min(const double& delEpsv) const;

  //--------------------------------------------------------------------------------------
  // Compute the function X2(kappa_{n+1})
  //  where
  //       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //--------------------------------------------------------------------------------------
  double
  computeX2(const double& kappa_old,
            const double& kappa_new,
            const double& G,
            const double& H,
            const double& delEpsv) const;

  //--------------------------------------------------------------------------------------
  // Compute the function dX2/dkappa(kappa_{n+1})
  //  where
  //       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
  //       Delta epsv = 0
  //--------------------------------------------------------------------------------------
  double
  computeDerivX2dkappa(const double& kappa_old,
                       const double& kappa_new,
                       const double& delEpsv) const;

  //--------------------------------------------------------------------------------------
  // Compute the constant B
  //  where
  //        B = 3 B0 [exp(p3+p4) - 1]
  //--------------------------------------------------------------------------------------
  double
  computeB() const;

  //--------------------------------------------------------------------------------------
  // Compute the function G(epsv)
  //  where
  //        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
  //        B = 3 B0 [exp(p3+p4) - 1]
  //        g34 = exp(p3+p4+epsv)
  //--------------------------------------------------------------------------------------
  double
  computeG(const double& epsv, const double& B) const;

  //--------------------------------------------------------------------------------------
  // Compute the function H(epsv)
  //  where
  //        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
  //        B = 3 B0 [exp(p3+p4) - 1]
  //        h3 = exp(p3+epsv)
  //--------------------------------------------------------------------------------------
  double
  computeH(const double& epsv, const double& B) const;

  //--------------------------------------------------------------------------------------
  // Compute the function F1(kappa, epsv) = f1(kappa) - G(epsv) + H(epsv)
  //  where f1(kappa) = 1/(p1 p3) exp(-p1 kappa - p0)
  //        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
  //        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
  //        B = 3 B0 [exp(p3+p4) - 1]
  //        g34 = exp(p3+p4+epsv)
  //        h3 = exp(p3+epsv)
  //--------------------------------------------------------------------------------------
  double
  computeF1(const double& kappa, const double& G, const double& H) const;

  //--------------------------------------------------------------------------------------
  // Compute the function dF1/dkappa(kappa, epsv) = df1/dkappa(kappa)
  //  where f1(kappa) = 1/(p1 p3) exp(-p1 kappa - p0)
  //--------------------------------------------------------------------------------------
  double
  computeDerivF1dkappa(const double& kappa) const;

  //--------------------------------------------------------------------------------------
  // Compute the function F2(kappa, epsv) = f2(kappa) - G(epsv) + H(epsv)
  //  where f2(kappa) = 1/(p1 p3) [kappa/p0]^(1-p0p1p3)
  //        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
  //        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
  //        B = 3 B0 [exp(p3+p4) - 1]
  //        g34 = exp(p3+p4+epsv)
  //        h3 = exp(p3+epsv)
  //--------------------------------------------------------------------------------------
  double
  computeF2(const double& kappa, const double& G, const double& H) const;

  //--------------------------------------------------------------------------------------
  // Compute the function dF2/dkappa(kappa, epsv) = df2/dkappa(kappa)
  //  where f2(kappa) = 1/(p1 p3) [kappa/p0]^(1-p0p1p3)
  //--------------------------------------------------------------------------------------
  double
  computeDerivF2dkappa(const double& kappa) const;
};

} // End namespace Uintah

#endif // __SOIL_MODEL_BRANNON_KAPPA_INT_VAR_MODEL_H__
