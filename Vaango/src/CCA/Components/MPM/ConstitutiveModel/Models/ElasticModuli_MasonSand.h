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

#ifndef ___ELASTIC_MODULI_MASON_SAND_MODEL_H__
#define ___ELASTIC_MODULI_MASON_SAND_MODEL_H__


#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Arenisca3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  /*! \class ElasticModuli_MasonSand
   *  \brief The elasticity from Micahel Homel's version of Arenisca3
   *  \author Biswajit Banerjee, 
   *
  */
  class ElasticModuli_MasonSand : public ElasticModuliModel {

  private:
    
    /* Tangent bulk modulus parameters */
    struct BulkModulusParameters {
      double b0;
      double b1;
      double b2;
      double Kmax;
      double alpha0;
      double alpha1;
      double alpha2;
      double alpha3;
      double alpha4;
    };

    /* Tangent shear modulus parameters */
    struct ShearModulusParameters {
      double G0;
      double G1;
      double G2;
      double G3;
      double G4;
    };

    /* Fluid bulk modulus parameters */
    struct FluidModulusParameters {
      double Kw;    // Bulk modulus of water
      double p0;    // Initial water pressure
      double gamma; // Air gamma = Cp/Cv
      double pRef;  // Reference air pressure (101325 Pa)
    };

    BulkModulusParameters d_bulk;
    ShearModulusParameters d_shear;
    FluidModulusParameters d_fluid;

    void checkInputParameters();

    ElasticModuli_MasonSand& operator=(const ElasticModuli_MasonSand &smm);

  public:
         
    /*! Construct a constant elasticity model. */
    ElasticModuli_MasonSand(Uintah::ProblemSpecP& ps);

    /*! Construct a copy of constant elasticity model. */
    ElasticModuli_MasonSand(const ElasticModuli_MasonSand* smm);

    /*! Destructor of constant elasticity model.   */
    virtual ~ElasticModuli_MasonSand();
         
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    /*! Compute the elasticity */
    ElasticModuli getInitialElasticModuli() const;
    ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const;
    ElasticModuli getElasticModuliLowerBound() const;
    ElasticModuli getElasticModuliUpperBound() const;

  };
} // End namespace Uintah
      
#endif  // __ELASTIC_MODULI_MASON_SAND_MODEL_H__

