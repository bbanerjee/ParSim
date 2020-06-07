/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#ifndef __BB_NO_KINEMATIC_HARDENING_MODEL_H__
#define __BB_NO_KINEMATIC_HARDENING_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardeningModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class KinematicHardening_None
  \brief Default kinematic hardening model - no kinematic hardening
  \author Biswajit Banerjee,
  Department of Mechanical Engineering,
  University of Utah

*/
/////////////////////////////////////////////////////////////////////////////

class KinematicHardening_None : public KinematicHardeningModel
{

private:
  // Prevent copying of this class
  // copy constructor
  // KinematicHardening_None(const KinematicHardening_None &cm);
  KinematicHardening_None& operator=(const KinematicHardening_None& cm);

public:
  // constructors
  KinematicHardening_None();
  KinematicHardening_None(Uintah::ProblemSpecP& ps);
  KinematicHardening_None(const KinematicHardening_None* cm);

  // destructor
  ~KinematicHardening_None() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["None"] = 0.0;
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

#endif // __BB_NO_KINEMATIC_HARDENING_MODEL_H__
