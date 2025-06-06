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

#ifndef ___ARENISCA_POISSON_RATIO_MODEL_H__
#define ___ARENISCA_POISSON_RATIO_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arenisca3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*! \class ElasticModuli_Arenisca
 *  \brief The elasticity from Micahel Homel's version of Arenisca3
 *  \author Biswajit Banerjee,
 *
*/
class ElasticModuli_Arenisca : public ElasticModuliModel
{

private:
  /* Tangent bulk modulus parameters */
  struct BulkModulusParameters
  {
    double B0;
    double B01;
    double B1;
    double B2;
    double B3;
    double B4;
  };

  /* Tangent shear modulus parameters */
  struct ShearModulusParameters
  {
    double G0;
    double G1;
    double G2;
    double G3;
    double G4;
  };

  BulkModulusParameters d_bulk;
  ShearModulusParameters d_shear;

  ElasticModuli_Arenisca& operator=(const ElasticModuli_Arenisca& smm);

public:
  /*! Construct a constant elasticity model. */
  ElasticModuli_Arenisca(Uintah::ProblemSpecP& ps);

  /*! Construct a copy of constant elasticity model. */
  ElasticModuli_Arenisca(const ElasticModuli_Arenisca* smm);

  /*! Destructor of constant elasticity model.   */
  ~ElasticModuli_Arenisca() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["B0"] = d_bulk.B0;
    params["B01"] = d_bulk.B01;
    params["B1"] = d_bulk.B1;
    params["B2"] = d_bulk.B2;
    params["B3"] = d_bulk.B3;
    params["B4"] = d_bulk.B4;
    params["G0"] = d_shear.G0;
    params["G1"] = d_shear.G1;
    params["G2"] = d_shear.G2;
    params["G3"] = d_shear.G3;
    params["G4"] = d_shear.G4;
    return params;
  }

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const override;
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const override;
  ElasticModuli getElasticModuliLowerBound() const override;
  ElasticModuli getElasticModuliUpperBound() const override;

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> computeDModuliDIntVar(const ModelStateBase* state) const override;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const override;
};
} // End namespace Vaango

#endif // __ARENISCA_POISSON_RATIO_MODEL_H__
