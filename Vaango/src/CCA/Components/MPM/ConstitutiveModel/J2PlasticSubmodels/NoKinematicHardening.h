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

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __NO_KINEMATIC_HARDENING_MODEL_H__
#define __NO_KINEMATIC_HARDENING_MODEL_H__

#include "KinematicHardeningModel.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class NoKinematicHardening
  \brief Default kinematic hardening model - no kinematic hardening
  \author Biswajit Banerjee,
  Department of Mechanical Engineering,
  University of Utah

*/
/////////////////////////////////////////////////////////////////////////////

class NoKinematicHardening : public KinematicHardeningModel
{

private:
  // Prevent copying of this class
  // copy constructor
  // NoKinematicHardening(const NoKinematicHardening &cm);
  NoKinematicHardening& operator=(const NoKinematicHardening& cm);

public:
  // constructors
  NoKinematicHardening();
  NoKinematicHardening(ProblemSpecP& ps);
  NoKinematicHardening(const NoKinematicHardening* cm);

  // destructor
  ~NoKinematicHardening() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  //////////
  /*! \brief Calculate the back stress */
  //////////
  void computeBackStress(const ModelStateBase* state, const double& delT,
                         const particleIndex idx, const double& delLambda,
                         const Matrix3& df_dsigma_new,
                         const Matrix3& backStress_old,
                         Matrix3& backStress_new) override;

  void eval_h_beta(const Matrix3& df_dsigma, const ModelStateBase* state,
                   Matrix3& h_beta) override;
};

} // End namespace Uintah

#endif // __NO_KINEMATIC_HARDENING_MODEL_H__
