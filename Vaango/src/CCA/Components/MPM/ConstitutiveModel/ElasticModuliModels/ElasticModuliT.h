/*
 * The MIT License
 *
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

#ifndef __ELASTICITY_TEMPLATED_MODEL_H__
#define __ELASTICITY_TEMPLATED_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateT.h>
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

/*! \class ElasticModuliT
 *  \brief A generic wrapper for various elasticity models
 *  \author Biswajit Banerjee,
 *
 * Provides a CRTP base class for various isotropic elasticity models
*/

template <typename DerivedT, typename StateT, typename PressureT, typename ShearT>
class ElasticModuliT
{

public:

  ~ElasticModuliT() = default;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps)
  {
    derived()->l_outputProblemSpec(ps);
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the model parameters
   */
  /////////////////////////////////////////////////////////////////////////
  ParameterDict
  getParameters() const
  {
    return derived()->l_getParameters();
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the elastic moduli
  */
  /////////////////////////////////////////////////////////////////////////
  ElasticModuli
  getInitialElasticModuli() const
  {
    return derived()->l_getInitialElasticModuli();
  }

  ElasticModuli
  getCurrentElasticModuli(const StateT* state) const
  {
    return derived()->l_getCurrentElasticModuli(state);
  }

  ElasticModuli
  getElasticModuliLowerBound() const
  {
    return derived()->l_getElasticModuliLowerBound();
  }

  ElasticModuli
  getElasticModuliUpperBound() const
  {
    return derived()->l_getElasticModuliUpperBound();
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the elastic moduli and their derivatives with respect to a single
           plastic internal variable
  */
  /////////////////////////////////////////////////////////////////////////
  std::pair<ElasticModuli, ElasticModuli>
  getElasticModuliAndDerivatives(const StateT* state) const
  {
    return derived()->l_getElasticModuliAndDerivatives(state);
  }

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> 
  computeDModuliDIntVar(const ModelStateBase* state) const
  {
    return derived()->l_computeDModuliDIntVar(state);
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief For partially saturated materials, get the drained and
           partially saturated moduli
  */
  /////////////////////////////////////////////////////////////////////////
  void
  computeDrainedModuli(const double& I1_bar,
                       const double& ev_p_bar,
                       double& KK,
                       double& GG)
  {
    derived()->l_computeDrainedModuli(I1_bar, ev_p_bar, KK, GG);
  }

  void
  computePartialSaturatedModuli(const double& I1_eff_bar,
                                const double& pw_bar,
                                const double& ev_p_bar,
                                const double& phi,
                                const double& S_w,
                                double& KK,
                                double& GG)
  {
    derived()->l_computePartialSaturatedModuli(I1_eff_bar, pw_bar, ev_p_bar, 
                                             phi, S_w, KK, GG);
  }

private:

  std::unique_ptr<PressureT> d_eos;
  std::unqiue_ptr<ShearT> d_shear;

  ElasticModuliT()
  {
  }

  DerivedT*
  derived()
  {
    return static_cast<DerivedT*>(this);
  }

  const DerivedT*
  derived() const
  {
    return static_cast<const DerivedT*>(this);
  }

  friend DerivedT;
};
} // End namespace Uintah

#endif // __ELASTICITY_TEMPLATED_MODEL_H__
