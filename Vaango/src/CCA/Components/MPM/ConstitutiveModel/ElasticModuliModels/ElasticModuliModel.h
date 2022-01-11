/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#ifndef __ELASTICITY_MODEL_H__
#define __ELASTICITY_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/TensorUtils.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

using ParameterDict = std::map<std::string, double>;

/*! \class ElasticModuli
 *  \brief A struct containing the instantaneous bulk and shear modulus
 *  \author Biswajit Banerjee,
*/
struct ElasticModuli
{
  double bulkModulus;
  double shearModulus;

  ElasticModuli(const double& bulk, const double& shear)
    : bulkModulus(bulk)
    , shearModulus(shear)
  {
  }
};

/*! \class ElasticModuliModel
 *  \brief A generic wrapper for various elasticity models
 *  \author Biswajit Banerjee,
 *
 * Provides an abstract base class for various isotropic elasticity models
*/

class ElasticModuliModel
{

public:
  //! Construct a elasticity model and initialize value
  /*! This is an abstract base class. */
  ElasticModuliModel();

  //! Destructor of elasticity model.
  /*! Virtual to ensure correct behavior */
  virtual ~ElasticModuliModel();

  virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the model parameters
   */
  /////////////////////////////////////////////////////////////////////////
  virtual ParameterDict getParameters() const = 0;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the elastic moduli
  */
  /////////////////////////////////////////////////////////////////////////
  virtual ElasticModuli getInitialElasticModuli() const = 0;
  virtual 
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const = 0; 

  virtual ElasticModuli getElasticModuliLowerBound() const = 0;
  virtual ElasticModuli getElasticModuliUpperBound() const = 0;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the elastic moduli and their derivatives with respect to a single
           plastic internal variable
  */
  /////////////////////////////////////////////////////////////////////////
  virtual std::pair<ElasticModuli, ElasticModuli>
    getElasticModuliAndDerivatives(const ModelStateBase* state_input) const 
  {
    return std::make_pair(ElasticModuli(0, 0), ElasticModuli(0, 0));
  }

  /*! Compute derivatives of moduli with respect to internal variables */
  virtual std::vector<ElasticModuli> 
  computeDModuliDIntVar(const ModelStateBase* state) const = 0;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  virtual std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const = 0;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief For partially saturated materials, get the drained and
           partially saturated moduli
  */
  /////////////////////////////////////////////////////////////////////////
  virtual void computeDrainedModuli(const double& I1_bar,
                                    const double& ev_p_bar, double& KK,
                                    double& GG) const {};
  virtual void computePartialSaturatedModuli(
    const double& I1_eff_bar, const double& pw_bar, const double& ev_p_bar,
    const double& phi, const double& S_w, double& KK, double& GG) const {};

  /*! Tangent modulus */
  Tensor::Matrix6Mandel
  computeElasticTangentModulus(const ModelStateBase* state) const
  {
    ElasticModuli moduli = getCurrentElasticModuli(state);
    double K = moduli.bulkModulus;
    double G = moduli.shearModulus;
    double K43G = K + 4.0 * G / 3.0;
    double K23G = K - 2.0 * G / 3.0;

    Tensor::Matrix6Mandel C_e = Tensor::Matrix6Mandel::Zero();
    C_e(0, 0) = K43G;
    C_e(0, 1) = K23G;
    C_e(0, 2) = K23G;
    C_e(1, 0) = K23G;
    C_e(1, 1) = K43G;
    C_e(1, 2) = K23G;
    C_e(2, 0) = K23G;
    C_e(2, 1) = K23G;
    C_e(2, 2) = K43G;
    C_e(3, 3) = G;
    C_e(4, 4) = G;
    C_e(5, 5) = G;
    return C_e;
  }

};
} // End namespace Vaango

#endif // __ELASTICITY_MODEL_H__
