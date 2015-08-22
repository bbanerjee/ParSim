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


#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_MasonSand.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Uintah;
using namespace Vaango;
         
// Construct a default elasticity model.  
/* From Arenisca3:
 * If the user has specified a nonzero G1 and G2, these are used to define a pressure
 * dependent poisson ratio, which is used to adjust the shear modulus along with the
 * bulk modulus.  The high pressure limit has nu=G1+G2;
 * if ((d_cm.G1!=0.0)&&(d_cm.G2!=0.0)){
 *     // High pressure bulk modulus:
 *     double nu = d_cm.G1+d_cm.G2;
 *     shear = 1.5*bulk*(1.0-2.0*nu)/(1.0+nu);
 * } 
 */

ElasticModuli_MasonSand::ElasticModuli_MasonSand(Uintah::ProblemSpecP& ps)
{
  ps->require("b0",     d_bulk.b0);      // Tangent Elastic Bulk Modulus Parameter
  ps->require("b1",     d_bulk.b1);      // Tangent Elastic Bulk Modulus Parameter
  ps->require("b2",     d_bulk.b2);      // Tangent Elastic Bulk Modulus Parameter
  ps->require("K_max",  d_bulk.Kmax);    // Tangent Elastic Bulk Modulus Parameter
  ps->require("alpha0", d_bulk.alpha0);  // Tangent Elastic Bulk Modulus Parameter
  ps->require("alpha1", d_bulk.alpha1);  // Tangent Elastic Bulk Modulus Parameter
  ps->require("alpha2", d_bulk.alpha2);  // Tangent Elastic Bulk Modulus Parameter
  ps->require("alpha3", d_bulk.alpha3);  // Tangent Elastic Bulk Modulus Parameter
  ps->require("alpha4", d_bulk.alpha4);  // Tangent Elastic Bulk Modulus Parameter

  ps->require("G0", d_shear.G0);         // Tangent Elastic Shear Modulus Parameter
  ps->require("G1", d_shear.G1);         // Tangent Elastic Shear Modulus Parameter
  ps->require("G2", d_shear.G2);         // Tangent Elastic Shear Modulus Parameter
  ps->require("G3", d_shear.G3);         // Tangent Elastic Shear Modulus Parameter
  ps->require("G4", d_shear.G4);         // Tangent Elastic Shear Modulus Parameter

  ps->require("K_w",   d_fluid.Kw);      // Bulk modulus of water
  ps->require("p0",    d_fluid.p0);      // Initial water pressure
  ps->require("gamma", d_fluid.gamma);   // Air gamma = Cp/Cv
  ps->require("p_ref", d_fluid.pRef);    // Reference air pressure (101325 Pa)

  checkInputParameters();
}

// Construct a copy of a elasticity model.  
ElasticModuli_MasonSand::ElasticModuli_MasonSand(const ElasticModuli_MasonSand* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

// Destructor of elasticity model.  
ElasticModuli_MasonSand::~ElasticModuli_MasonSand()
{
}

void 
ElasticModuli_MasonSand::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "mason_sand");

  elasticModuli_ps->appendElement("b0",d_bulk.b0);
  elasticModuli_ps->appendElement("b1",d_bulk.b1);
  elasticModuli_ps->appendElement("b2",d_bulk.b2);
  elasticModuli_ps->appendElement("K_max",d_bulk.Kmax);
  elasticModuli_ps->appendElement("alpha0",d_bulk.alpha0);
  elasticModuli_ps->appendElement("alpha1",d_bulk.alpha1);
  elasticModuli_ps->appendElement("alpha2",d_bulk.alpha2);
  elasticModuli_ps->appendElement("alpha3",d_bulk.alpha3);
  elasticModuli_ps->appendElement("alpha4",d_bulk.alpha4);

  elasticModuli_ps->appendElement("G0",d_shear.G0);
  elasticModuli_ps->appendElement("G1",d_shear.G1);  // Low pressure Poisson ratio
  elasticModuli_ps->appendElement("G2",d_shear.G2);  // Pressure-dependent Poisson ratio term
  elasticModuli_ps->appendElement("G3",d_shear.G3);  // Not used
  elasticModuli_ps->appendElement("G4",d_shear.G4);  // Not used

  elasticModuli_ps->appendElement("K_w",   d_fluid.Kw);    // Bulk modulus of water
  elasticModuli_ps->appendElement("p0",    d_fluid.p0);    // Initial water pressure
  elasticModuli_ps->appendElement("gamma", d_fluid.gamma); // Air gamma = Cp/Cv
  elasticModuli_ps->appendElement("p_ref", d_fluid.pRef);  // Refernec air pressure (101325 Pa)
}
         
// Compute the elasticity
ElasticModuli 
ElasticModuli_MasonSand::getInitialElasticModuli() const
{
  return ElasticModuli(d_bulk.Kmax*d_bulk.b0, d_shear.G0);
}

ElasticModuli 
ElasticModuli_MasonSand::getElasticModuliUpperBound() const
{
  return ElasticModuli(d_bulk.Kmax*(d_bulk.b0+d_bulk.b1), d_shear.G0);
}

ElasticModuli 
ElasticModuli_MasonSand::getElasticModuliLowerBound() const
{
  return ElasticModuli(d_bulk.Kmax*d_bulk.b0, d_shear.G0);
}

ElasticModuli 
ElasticModuli_MasonSand::getCurrentElasticModuli(const ModelStateBase* state_input) const
{
  const ModelState_Arenisca3* state = dynamic_cast<const ModelState_Arenisca3*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Make sure the quantities are positive in compression
  double I1_bar = -state->I1;
  double ev_p_bar = -(state->plasticStrainTensor).Trace();  

  double KK = d_bulk.b0*d_bulk.Kmax;
  double GG = d_shear.G0;
  if (I1_bar > 0.0) {
    
    double I1_sat = d_bulk.alpha0;
    if (ev_p_bar > 0) {
      I1_sat += d_bulk.alpha1/(d_bulk.alpha2 + d_bulk.alpha3*std::exp(-d_bulk.alpha4*ev_p_bar));
    } else {
      I1_sat += d_bulk.alpha1/(d_bulk.alpha2 + d_bulk.alpha3);
    }

    KK += d_bulk.Kmax*d_bulk.b1*std::pow(I1_bar/I1_sat, d_bulk.b2);
    
    double nu = d_shear.G1 + d_shear.G2*exp(-I1_bar/I1_sat);
    GG = (nu > 0.0) ? 1.5*KK*(1.0-2.0*nu)/(1.0+nu) : GG;
  } 

  return ElasticModuli(KK, GG);
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
ElasticModuli_MasonSand::checkInputParameters()
{
  std::ostringstream warn;

  if (d_bulk.b0 <= 0.0) {
    warn << "b0 must be positive. b0 = " << d_bulk.b0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk.b1 < 0.0) {
    warn << "b1 must be nonnegative. b1 = " << d_bulk.b1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk.b2 < 0.0) {
    warn << "b2 must be nonnegative. b2 = " << d_bulk.b2 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_fluid.Kw < 0.0) {
    warn << "Fluid Kw must be >=0. Kw = " << d_fluid.Kw << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_fluid.p0 < 0.0) {
    warn << "Negative initial fluid pressure p0 not supported. p0 = " << d_fluid.p0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_fluid.Kw < 0.0 && (d_bulk.b0 == 0.0 || d_bulk.b1 == 0.0)) {
    warn << "Solid b0 and b1 must be positive to use fluid model." << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  
}

