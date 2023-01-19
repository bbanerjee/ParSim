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

#include <CCA/Components/MPM/ParticleCreator/FractureParticleCreator.h>
#include <Core/Grid/Box.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/CrackBC.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>

using namespace Uintah;
using std::vector;

FractureParticleCreator::FractureParticleCreator(MPMMaterial* matl,
                                                 MPMFlags* flags)

  :  ParticleCreator(matl,flags)
{
  registerPermanentParticleState(matl);
}

FractureParticleCreator::~FractureParticleCreator()
{
}

void
FractureParticleCreator::registerPermanentParticleState(MPMMaterial* /*matl*/)

{
  //particle_state.push_back(lb->pX0Label);
  //particle_state_preReloc.push_back(lb->pX0Label_preReloc);

}


void 
FractureParticleCreator::applyForceBC(const Vector& dxpp, 
                                      const Point& pp,
                                      const double& pMass, 
                                      Vector& pExtForce)
{
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();
        
    //cerr << " BC Type = " << bcType << endl;
    if (bcType == "Force") {
      ForceBC* fbc = dynamic_cast<ForceBC*>(bc.get());

      Box fbcBox(fbc->getLowerRange(), fbc->getUpperRange());

      //cerr << "BC Box = " << bcBox << " Point = " << pp << endl;
      if(fbcBox.contains(pp)) {
        pExtForce = fbc->getForceDensity() * pMass;
        //cerr << "External Force on Particle = " << pExtForce 
        //     << " Force Density = " << fbc->getForceDensity() 
        //     << " Particle Mass = " << pMass << endl;
      }
    } 
  }
}
