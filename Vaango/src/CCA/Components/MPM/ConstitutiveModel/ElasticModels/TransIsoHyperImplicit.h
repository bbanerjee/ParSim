/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Reseach Limited, New Zealand
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

#ifndef __Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__
#define __Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyper.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <Core/Math/Matrix3.h>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace Uintah {

using Matrix66 = Eigen::Matrix<double, 6, 6>;

class TransIsoHyperImplicit : public TransIsoHyper, public ImplicitCM
{
public:

  TransIsoHyperImplicit(ProblemSpecP& ps, MPMFlags* flag);
  TransIsoHyperImplicit(const TransIsoHyperImplicit* cm);
  TransIsoHyperImplicit& operator=(const TransIsoHyperImplicit& cm) = delete;
  ~TransIsoHyperImplicit() override;
  
  std::unique_ptr<ConstitutiveModel> clone() override;

  ModelType modelType() const override { return ModelType::TOTAL_FORM; }

  /* initialization */
  void initializeCMData(const Patch* patch,
                        const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  /* For newton iterations */
  void addComputesAndRequires(Task* task,
                              const MPMMaterial* matl,
                              const PatchSet* patches,
                              const bool recursion,
                              const bool SchedParent) const override;
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   Solver* solver,
                                   const bool recursion) override;

  /* For end of iterations */
  void addComputesAndRequires(Task* task,
                              const MPMMaterial* matl,
                              const PatchSet* patches) const override;
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw) override;

private:

  Matrix66 computeTangentStiffness(const TransIsoHyperState& ss,
                                   double p,
                                   double pFailed) const;
  Matrix66 computeVolTangentStiffness(const TransIsoHyperState& ss, double p) const;
  Matrix66 computeDevTangentStiffness(const TransIsoHyperState& ss) const;
  Matrix66 computeFiberTangentStiffness(const TransIsoHyperState& ss) const;

};
} // End namespace Uintah

#endif // __Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__
