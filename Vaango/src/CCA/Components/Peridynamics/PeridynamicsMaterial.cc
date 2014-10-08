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

//  PeridynamicsMaterial.cc

#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModelFactory.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModelFactory.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreatorFactory.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/SimulationState.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/GeometryPiece/NullGeometryPiece.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/PerPatch.h>

#include <Core/Util/DebugStream.h>

#include <iostream>
#include <string>
#include <list>

using namespace Vaango;

using Uintah::DebugStream;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDMatDoing:+,PDMatDebug:+".....
//  bash     : export SCI_DEBUG="PDMatDoing:+,PDMatDebug:+" )
//  default is OFF
static DebugStream cout_doing("PDMatDoing", false);
static DebugStream cout_dbg("PDMatDebug", false);

PeridynamicsMaterial::PeridynamicsMaterial(Uintah::ProblemSpecP& ps, 
                                           Uintah::SimulationStateP& ss,
                                           PeridynamicsFlags* flags)
  : Uintah::Material(ps), d_materialModel(0), d_particle_creator(0)
{
  d_varLabel = scinew PeridynamicsLabel();

  // The standard set of initializations needed
  standardInitialization(ps, flags);
  
  d_materialModel->setSharedState(ss.get_rep());
  //d_damageModel->setSharedState(ss.get_rep());  // TODO: Not sure about this

  // Check to see which ParticleCreator object we need
  d_particle_creator = ParticleCreatorFactory::create(ps, this, flags);
}

void
PeridynamicsMaterial::standardInitialization(Uintah::ProblemSpecP& ps, 
                                             PeridynamicsFlags* flags)
{
  ps->require("density", d_density);

  d_materialModel = PeridynamicsMaterialModelFactory::create(ps,flags);
  if(!d_materialModel){
    std::ostringstream desc;
    desc << "An error occured in the PeridynamicsMaterialModelFactory that has \n" 
         << " slipped through the existing bullet proofing. " << std::endl;
    throw Uintah::ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_damageModel = PeridynamicsDamageModelFactory::create(ps, d_varLabel, flags);
  if(!d_damageModel){
    std::ostringstream desc;
    desc << "An error occured in the PeridynamicsDamageModelFactory that has \n" 
         << " slipped through the existing bullet proofing. " << std::endl;
    throw Uintah::ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  std::list<Uintah::GeometryObject::DataItem> geom_obj_data;
  geom_obj_data.push_back(Uintah::GeometryObject::DataItem("res", Uintah::GeometryObject::IntVector));
  geom_obj_data.push_back(Uintah::GeometryObject::DataItem("temperature", Uintah::GeometryObject::Double));
  geom_obj_data.push_back(Uintah::GeometryObject::DataItem("velocity", Uintah::GeometryObject::Vector));

  for (Uintah::ProblemSpecP geom_obj_ps = ps->findBlock("geom_object");
       geom_obj_ps != 0; 
       geom_obj_ps = geom_obj_ps->findNextBlock("geom_object") ) {

    std::vector<Uintah::GeometryPieceP> pieces;
    Uintah::GeometryPieceFactory::create(geom_obj_ps, pieces);

    Uintah::GeometryPieceP mainpiece;
    if(pieces.size() == 0){
      throw Uintah::ParameterNotFound("No piece specified in geom_object", __FILE__, __LINE__);
    } else if(pieces.size() > 1){
      mainpiece = scinew Uintah::UnionGeometryPiece(pieces);
    } else {
      mainpiece = pieces[0];
    }

    //    piece_num++;
    d_geom_objs.push_back(scinew Uintah::GeometryObject(mainpiece, geom_obj_ps, geom_obj_data));
  }
}

// Default constructor
PeridynamicsMaterial::PeridynamicsMaterial() : 
  d_materialModel(0), d_particle_creator(0)
{
  d_varLabel = scinew PeridynamicsLabel();
  d_damageModel = 0;
}

PeridynamicsMaterial::~PeridynamicsMaterial()
{
  delete d_varLabel;
  delete d_materialModel;
  delete d_particle_creator;
  delete d_damageModel;

  for (int i = 0; i<(int)d_geom_objs.size(); i++) {
    delete d_geom_objs[i];
  }
}

/*
*/
void PeridynamicsMaterial::registerParticleState(Uintah::SimulationState* sharedState)
{
  sharedState->d_particleState.push_back(d_particle_creator->returnParticleState());
  sharedState->d_particleState_preReloc.push_back(d_particle_creator->returnParticleStatePreReloc());
}

Uintah::ProblemSpecP 
PeridynamicsMaterial::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP peridynamics_ps = Uintah::Material::outputProblemSpec(ps);
  d_materialModel->outputProblemSpec(peridynamics_ps);
  d_damageModel->outputProblemSpec(peridynamics_ps);
  for (std::vector<Uintah::GeometryObject*>::const_iterator it = d_geom_objs.begin();
       it != d_geom_objs.end(); it++) {
    (*it)->outputProblemSpec(peridynamics_ps);
  }

  return peridynamics_ps;
}

PeridynamicsMaterialModel* 
PeridynamicsMaterial::getMaterialModel() const
{
  return d_materialModel;
}

PeridynamicsDamageModel* 
PeridynamicsMaterial::getDamageModel() const
{
  return d_damageModel;
}

Uintah::particleIndex 
PeridynamicsMaterial::countParticles(const Uintah::Patch* patch)
{
  return d_particle_creator->countParticles(patch, d_geom_objs);
}

void 
PeridynamicsMaterial::createParticles(Uintah::particleIndex numParticles,
                                      Uintah::CCVariable<short int>& cellNAPID,
                                      const Uintah::Patch* patch,
                                      Uintah::DataWarehouse* new_dw)
{
  d_particle_creator->createParticles(this,numParticles,cellNAPID,
                                      patch,new_dw,d_geom_objs);
}

ParticleCreator* 
PeridynamicsMaterial::getParticleCreator()
{
  return  d_particle_creator;
}


