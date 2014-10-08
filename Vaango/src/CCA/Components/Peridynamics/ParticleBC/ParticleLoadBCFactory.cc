/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCFactory.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Vaango;

std::vector<ParticleLoadBCBase*> ParticleLoadBCFactory::particleLoadBCs;

void 
ParticleLoadBCFactory::create(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP test = ps->findBlock("ParticleBC");
  if (test){

    Uintah::ProblemSpecP current_ps = ps->findBlock("ParticleBC")->findBlock("Load");


    for(Uintah::ProblemSpecP child = current_ps->findBlock("force"); child != 0;
        child = child->findNextBlock("force") ) {
       particleLoadBCs.push_back(scinew ParticleForceBC(child));
    }

    for(Uintah::ProblemSpecP child = current_ps->findBlock("normal_force"); child != 0;
        child = child->findNextBlock("normal_force") ) {
       particleLoadBCs.push_back(scinew ParticleNormalForceBC(child));
    }

    for(Uintah::ProblemSpecP child = current_ps->findBlock("pressure"); child != 0;
        child = child->findNextBlock("pressure") ) {
       particleLoadBCs.push_back(scinew ParticlePressureBC(child));
    }

  }
}

void ParticleLoadBCFactory::clean()
{
  for (int i = 0; i < static_cast<int>(particleLoadBCs.size()); i++)
    delete particleLoadBCs[i];
  particleLoadBCs.clear();
}
