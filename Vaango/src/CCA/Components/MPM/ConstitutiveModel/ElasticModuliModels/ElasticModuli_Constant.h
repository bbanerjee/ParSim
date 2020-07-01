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

#ifndef __BB_CONSTANT_ELASTICITY_MODEL_H__
#define __BB_CONSTANT_ELASTICITY_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*! \class ElasticModuli_Constant
 *  \brief The elasticity does not vary with density and temperature
 *  \author Biswajit Banerjee,
 *
*/
class ElasticModuli_Constant : public ElasticModuliModel
{

private:
  double d_bulk;
  double d_shear;

  ElasticModuli_Constant& operator=(const ElasticModuli_Constant& smm);

public:
  /*! Construct a constant elasticity model. */
  ElasticModuli_Constant(Uintah::ProblemSpecP& ps);

  /*! Construct a copy of constant elasticity model. */
  ElasticModuli_Constant(const ElasticModuli_Constant* smm);

  /*! Destructor of constant elasticity model.   */
  ~ElasticModuli_Constant() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["K"] = d_bulk;
    params["G"] = d_shear;
    return params;
  }

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const override;
  ElasticModuli getCurrentElasticModuli(const ModelStateBase*) const override;
  ElasticModuli getElasticModuliLowerBound() const override;
  ElasticModuli getElasticModuliUpperBound() const override;

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> computeDModuliDIntVar(const ModelStateBase* state) const override;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const override;
};
} // End namespace Vaango

#endif // __CONSTANT_ELASTICITY_MODEL_H__
