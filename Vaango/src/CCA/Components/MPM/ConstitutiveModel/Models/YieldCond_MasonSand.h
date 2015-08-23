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

#ifndef __PARTIALLy_SATURATED_ARENISCA3_YIELD_CONDITION_MODEL_H__
#define __PARTIALLy_SATURATED_ARENISCA3_YIELD_CONDITION_MODEL_H__


#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  /*! 
   \class  YieldCond_MasonSand
   \brief  The Partally saturated Arenisca3 yield condition
  */

  class YieldCond_MasonSand : public YieldCondition {
  
  friend class InternalVar_MasonSand;

  public:
    
    static const double sqrt_three;
    static const double one_sqrt_three;

  private:

    struct InputParameters {
      double PEAKI1;
      double FSLOPE;
      double STREN;
      double YSLOPE;
      double BETA;
    };

    struct ModelParameters {
      double a1;
      double a2;
      double a3;
      double a4;
      double beta;
      double CR;
    };

    struct CapParameters {
      double CR;
    };

    struct RateParameters {
      double T1;
      double T2;
    };

    InputParameters d_inputParam;
    ModelParameters d_modelParam;
    CapParameters   d_capParam;
    RateParameters  d_rateParam;

    void checkInputParameters();
    void computeModelParameters();

    // Prevent copying of this class
    // copy constructor
    //YieldCond_MasonSand(const YieldCond_MasonSand &);
    YieldCond_MasonSand& operator=(const YieldCond_MasonSand &);

  public:

    //! Constructor
    /*! Creates a YieldCond_MasonSand function object */
    YieldCond_MasonSand(Uintah::ProblemSpecP& ps);
    YieldCond_MasonSand(const YieldCond_MasonSand* cm);
         
    //! Destructor 
    ~YieldCond_MasonSand();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    /*! Get parameters */
    std::map<std::string, double> getParameters() const {
      std::map<std::string, double> params;
      params["PEAKI1"] = d_inputParam.PEAKI1;
      params["FSLOPE"] = d_inputParam.FSLOPE;
      params["STREN"] = d_inputParam.STREN;
      params["YSLOPE"] = d_inputParam.YSLOPE;
      params["BETA"] = d_inputParam.BETA;
      params["CR"] = d_capParam.CR;
      params["T1"] = d_rateParam.T1;
      params["T2"] = d_rateParam.T2;
      params["a1"] = d_modelParam.a1;
      params["a2"] = d_modelParam.a2;
      params["a3"] = d_modelParam.a3;
      params["a4"] = d_modelParam.a4;
      return params;
    }

    //--------------------------------------------------------------
    // Compute value of yield function
    //--------------------------------------------------------------
    double evalYieldCondition(const ModelStateBase* state);
    double evalYieldConditionMax(const ModelStateBase* state);

    //--------------------------------------------------------------
    // Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
    //--------------------------------------------------------------
    double computeVolStressDerivOfYieldFunction(const ModelStateBase* state);

    //--------------------------------------------------------------
    // Compute df/dq  where q = sqrt(3 J_2), J_2 = 2nd invariant deviatoric stress
    //--------------------------------------------------------------
    double computeDevStressDerivOfYieldFunction(const ModelStateBase* state);

    //--------------------------------------------------------------
    // Compute d/depse_v(df/dp)
    //--------------------------------------------------------------
    double computeVolStrainDerivOfDfDp(const ModelStateBase* state,
                                       const PressureModel* eos,
                                       const ShearModulusModel* shear,
                                       const InternalVariableModel* intvar);

    //--------------------------------------------------------------
    // Compute d/depse_s(df/dp)
    //--------------------------------------------------------------
    double computeDevStrainDerivOfDfDp(const ModelStateBase* state,
                                       const PressureModel* eos,
                                       const ShearModulusModel* shear,
                                       const InternalVariableModel* intvar);

    //--------------------------------------------------------------
    // Compute d/depse_v(df/dq)
    //--------------------------------------------------------------
    double computeVolStrainDerivOfDfDq(const ModelStateBase* state,
                                       const PressureModel* eos,
                                       const ShearModulusModel* shear,
                                       const InternalVariableModel* intvar);

    //--------------------------------------------------------------
    // Compute d/depse_s(df/dq)
    //--------------------------------------------------------------
    double computeDevStrainDerivOfDfDq(const ModelStateBase* state,
                                       const PressureModel* eos,
                                       const ShearModulusModel* shear,
                                       const InternalVariableModel* intvar);

    //--------------------------------------------------------------
    // Compute df/depse_v
    //--------------------------------------------------------------
    double computeVolStrainDerivOfYieldFunction(const ModelStateBase* state,
                                                const PressureModel* eos,
                                                const ShearModulusModel* shear,
                                                const InternalVariableModel* intvar);

    //--------------------------------------------------------------
    // Compute df/depse_s
    //--------------------------------------------------------------
    double computeDevStrainDerivOfYieldFunction(const ModelStateBase* state,
                                                const PressureModel* eos,
                                                const ShearModulusModel* shear,
                                                const InternalVariableModel* intvar);

    //================================================================================
    // Other options below.
    //================================================================================

    // Evaluate the yield function.
    double evalYieldCondition(const double p,
                              const double q,
                              const double dummy0,
                              const double dummy1,
                              double& dummy2);

    // Evaluate yield condition (s = deviatoric stress = sigDev
    //                           p = state->pressure
    //                           p_c = state->yieldStress)
    double evalYieldCondition(const Uintah::Matrix3& sigDev,
                              const ModelStateBase* state);

    /////////////////////////////////////////////////////////////////////////
    /*! 
      \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
      with respect to \f$\sigma_{ij}\f$.
    */
    /////////////////////////////////////////////////////////////////////////
    void evalDerivOfYieldFunction(const Uintah::Matrix3& stress,
                                  const double dummy1,
                                  const double dummy2,
                                  Uintah::Matrix3& derivative);

    /////////////////////////////////////////////////////////////////////////
    /*! 
      \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
      with respect to \f$s_{ij}\f$.

      This is for the associated flow rule with \f$s_{ij}\f$ being
      the deviatoric stress.
    */
    /////////////////////////////////////////////////////////////////////////
    void evalDevDerivOfYieldFunction(const Uintah::Matrix3& stress,
                                     const double dummy1,
                                     const double dummy2,
                                     Uintah::Matrix3& derivative);

    /*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
    void eval_df_dsigma(const Uintah::Matrix3& xi,
                        const ModelStateBase* state,
                        Uintah::Matrix3& df_dsigma);

    /*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$  
        where \f$s\f$ is deviatoric part of Cauchy stress and 
        \f$\beta\f$ is the backstress */
    void eval_df_dxi(const Uintah::Matrix3& xi,
                     const ModelStateBase* state,
                     Uintah::Matrix3& df_xi);

    /* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
    void eval_df_ds_df_dbeta(const Uintah::Matrix3& xi,
                             const ModelStateBase* state,
                             Uintah::Matrix3& df_ds,
                             Uintah::Matrix3& df_dbeta);

    /*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$)*/
    double eval_df_dep(const Uintah::Matrix3& xi,
                       const double& d_sigy_dep,
                       const ModelStateBase* state);

    /*! Derivative with respect to the porosity (\f$\epsilon^p \f$)*/
    double eval_df_dphi(const Uintah::Matrix3& xi,
                        const ModelStateBase* state);

    /*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
    double eval_h_alpha(const Uintah::Matrix3& xi,
                        const ModelStateBase* state);

    /*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
    double eval_h_phi(const Uintah::Matrix3& xi,
                      const double& factorA,
                      const ModelStateBase* state);

    /////////////////////////////////////////////////////////////////////////
    /*! 
      \brief Compute the elastic-plastic tangent modulus.
    */
    /////////////////////////////////////////////////////////////////////////
    void computeElasPlasTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                       const Uintah::Matrix3& sigma, 
                                       double sigY,
                                       double dsigYdep,
                                       double porosity,
                                       double voidNuclFac,
                                       Uintah::TangentModulusTensor& Cep);

    /////////////////////////////////////////////////////////////////////////
    /*!
      \brief Evaluate the factor \f$h_1\f$ for plastic strain

      \f[
      h_1 = \frac{\sigma : f_{\sigma}}{\sigma_Y}
      \f]
      
      \return factor 
    */
    /////////////////////////////////////////////////////////////////////////
    inline double computePlasticStrainFactor(double sigma_f_sigma,
                                             double sigma_Y);

    /////////////////////////////////////////////////////////////////////////
    /*!
      \brief Compute the continuum elasto-plastic tangent modulus
      assuming associated flow rule.

      \f[
      C_{ep} = C_{e} - \frac{(C_e:f_{\sigma})\otimes(f_{\sigma}:C_e)}
      {-f_q.h_q + f_{\sigma}:C_e:f_{\sigma}}
      \f]
      
      \return TangentModulusTensor \f$ C_{ep} \f$.
    */
    /////////////////////////////////////////////////////////////////////////
    void computeTangentModulus(const Uintah::TangentModulusTensor& Ce,
                               const Uintah::Matrix3& f_sigma, 
                               double f_q1, 
                               double h_q1,
                               Uintah::TangentModulusTensor& Cep);

  };

} // End namespace Uintah

#endif  // __PARTIALLy_SATURATED_ARENISCA3_YIELD_CONDITION_MODEL_H__ 
