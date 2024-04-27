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

#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>

#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/CrackBC.h>
#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/ArchesHeatFluxBC.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Uintah;

std::vector<MPMPhysicalBC_SP> MPMPhysicalBCFactory::mpmPhysicalBCs;

void 
MPMPhysicalBCFactory::create(const ProblemSpecP& ps, 
                             const GridP& grid, 
                             const MPMFlags* flags)
{
  ProblemSpecP test = ps->findBlock("PhysicalBC");
  if (test){

    ProblemSpecP current_ps = ps->findBlock("PhysicalBC")->findBlock("MPM");

    for(ProblemSpecP child = current_ps->findBlock("force"); child != 0;
        child = child->findNextBlock("force") ) {
      mpmPhysicalBCs.push_back(std::make_unique<ForceBC>(child));
    }

    for(ProblemSpecP child = current_ps->findBlock("moment"); child != 0;
        child = child->findNextBlock("moment") ) {
      mpmPhysicalBCs.push_back(std::make_unique<MomentBC>(child, grid, flags));
    }

    for(ProblemSpecP child = current_ps->findBlock("velocity"); child != 0;
        child = child->findNextBlock("velocity") ) {
      mpmPhysicalBCs.push_back(std::make_unique<VelocityBC>(child, grid, flags));
    }

    for(ProblemSpecP child = current_ps->findBlock("pressure"); child != 0;
        child = child->findNextBlock("pressure") ) {
      mpmPhysicalBCs.push_back(std::make_unique<PressureBC>(child, grid, flags));
    }

    for(ProblemSpecP child = current_ps->findBlock("crack"); child != 0;
        child = child->findNextBlock("crack") ) {
      mpmPhysicalBCs.push_back(std::make_unique<CrackBC>(child));
    }

    for(ProblemSpecP child = current_ps->findBlock("heat_flux"); child != 0;
        child = child->findNextBlock("heat_flux") ) {
      mpmPhysicalBCs.push_back(std::make_unique<HeatFluxBC>(child, grid));
    }

  }
}

void 
MPMPhysicalBCFactory::clean()
{
  mpmPhysicalBCs.clear();
}
