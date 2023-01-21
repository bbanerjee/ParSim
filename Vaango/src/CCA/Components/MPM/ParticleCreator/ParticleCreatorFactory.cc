/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ParticleCreator/ParticleCreatorFactory.h>

#include <CCA/Components/MPM/ParticleCreator/FractureParticleCreator.h>
#include <CCA/Components/MPM/ParticleCreator/ImplicitParticleCreator.h>
#include <CCA/Components/MPM/ParticleCreator/MembraneParticleCreator.h>
#include <CCA/Components/MPM/ParticleCreator/ShellParticleCreator.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <Core/ProblemSpec/ProblemSpec.h>

namespace Uintah {

std::unique_ptr<ParticleCreator>
ParticleCreatorFactory::create(ProblemSpecP& ps,
                               MPMMaterial* mat,
                               MPMFlags* flags)
{

  ProblemSpecP cm_ps = ps->findBlock("constitutive_model");
  string mat_type;
  cm_ps->getAttribute("type", mat_type);

  if (flags->d_integratorType == "implicit") {
    return std::make_unique<ImplicitParticleCreator>(mat, flags);
  }

  else if (flags->d_integratorType == "fracture") {
    return std::make_unique<FractureParticleCreator>(mat, flags);
  }

  else if (mat_type == "membrane") {
    return std::make_unique<MembraneParticleCreator>(mat, flags);
  }

  else if (mat_type == "shell_CNH") {
    return std::make_unique<ShellParticleCreator>(mat, flags);
  }

  else {
    return std::make_unique<ParticleCreator>(mat, flags);
  }
}

} // end namespace Uintah