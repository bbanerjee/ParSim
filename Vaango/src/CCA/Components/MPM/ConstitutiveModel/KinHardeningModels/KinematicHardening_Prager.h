/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __BB_PRAGER_KINEMATIC_HARDENING_MODEL_H__
#define __BB_PRAGER_KINEMATIC_HARDENING_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardeningModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class KinematicHardening_Prager
  \brief Default kinematic hardening model - no kinematic hardening
  \author Biswajit Banerjee,
  Department of Mechanical Engineering,
  University of Utah

  The kinematic hardening rule is given by
  \f[
  \dot{\beta} = \frac{2}{3}~H_prime~d_p
  \f]
  where \f$\beta\f$ is the back stress, \f$H_prime\f$ is the hardening
  modulus, and \f$ d_p\f$ is the plastic rate of deformation.

  For associative plasticity
  \f[
    d_p = dot{\lambda}~\frac{\partial f}{\partial sigma}
  \f]
  For von Mises plasticity with the yield condition of the form
  \f[
    f := \sqrt{\frac{3}{2} \xi:\xi} - \sigma_y \le 0
  \f]
  where \f$\xi = s - \beta\f$ and \f$s\f$ is the deviatoric part of the Cauchy
  stress, we have
  \f[
    \frac{\partial f}{\partial sigma} = \sqrt{\frac{3}{2}}\cfrac{\xi}{||\xi||}
     = \sqrt{\frac{3}{2}}~n ~;~~ ||n|| = 1
  \f]
  and
  \f[
    d_p = \sqrt{\frac{3}{2}}~\dot{\lambda}~n
  \f]
  Therefore, the evolution equation for beta can be written as
  \f[
  \dot{\beta} = \sqrt{\frac{2}{3}}~H_1~\dot{\lambda}~n
  \f]
  A backward Euler discretization leads to
  \f[
  \beta_{n+1} = \beta_n + \sqrt{\frac{2}{3}}~H_1~\Delta\lambda~n_{n+1}
  \f]
*/
/////////////////////////////////////////////////////////////////////////////

class KinematicHardening_Prager : public KinematicHardeningModel
{

protected:
  struct CMData
  {
    double beta;              // beta is a parameter between 0 and 1
                              // 0 == no kinematic hardening
    double hardening_modulus; // the kinematic hardening modulus
  };

private:
  CMData d_cm;

  // Prevent copying of this class
  // copy constructor
  // KinematicHardening_Prager(const KinematicHardening_Prager &cm);
  KinematicHardening_Prager& operator=(const KinematicHardening_Prager& cm);

public:
  // constructors
  KinematicHardening_Prager(Uintah::ProblemSpecP& ps);
  KinematicHardening_Prager(const KinematicHardening_Prager* cm);

  // destructor
  ~KinematicHardening_Prager() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["beta"] = d_cm.beta;
    params["H1"] = d_cm.hardening_modulus;
    return params;
  }

  //////////
  /*! \brief Calculate the back stress */
  //////////
  void computeBackStress(const ModelStateBase* state, const double& delT,
                         const Uintah::particleIndex idx,
                         const double& delLambda,
                         const Uintah::Matrix3& df_dsigma_new,
                         const Uintah::Matrix3& backStress_old,
                         Uintah::Matrix3& backStress_new) override;

  void eval_h_beta(const Uintah::Matrix3& df_dsigma,
                   const ModelStateBase* state,
                   Uintah::Matrix3& h_beta) override;
};

} // End namespace Uintah

#endif // __BB_PRAGER_KINEMATIC_HARDENING_MODEL_H__
