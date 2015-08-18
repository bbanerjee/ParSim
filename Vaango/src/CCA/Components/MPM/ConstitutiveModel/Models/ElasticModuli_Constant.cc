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


#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_Constant.h>

using namespace Uintah;
using namespace Vaango;
         
// Construct a default elasticity model.  
ElasticModuli_Constant::ElasticModuli_Constant(Uintah::ProblemSpecP& ps)
{
  ps->require("bulk_modulus", d_bulk);
  ps->require("shear_modulus", d_shear);
}

// Construct a copy of a elasticity model.  
ElasticModuli_Constant::ElasticModuli_Constant(const ElasticModuli_Constant* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

// Destructor of elasticity model.  
ElasticModuli_Constant::~ElasticModuli_Constant()
{
}

void ElasticModuli_Constant::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type","constant");
  elasticModuli_ps->appendElement("bulk_modulus", d_bulk);
  elasticModuli_ps->appendElement("shear_modulus", d_shear);
}
         
// Compute the elasticity
ElasticModuli
ElasticModuli_Constant::getInitialElasticModuli()
{
  return ElasticModuli(d_bulk, d_shear);
}

ElasticModuli
ElasticModuli_Constant::getCurrentElasticModuli(const ModelState* state) const
{
  return ElasticModuli(d_bulk, d_shear);
}

ElasticModuli 
ElasticModuli_Constant::getElasticModuliLowerBound() const
{
  return ElasticModuli(d_bulk, d_shear);
}

ElasticModuli 
ElasticModuli_Constant::getElasticModuliUpperBound() const
{
  return ElasticModuli(d_bulk, d_shear);
}

