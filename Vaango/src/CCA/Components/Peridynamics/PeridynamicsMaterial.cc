//  PeridynamicsMaterial.cc

#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModelFactory.h>
#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>
#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModelFactory.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreatorFactory.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>
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
#include <iostream>
#include <string>
#include <list>

using namespace Vaango;

PeridynamicsMaterial::PeridynamicsMaterial(Uintah::ProblemSpecP& ps, 
                                           Uintah::SimulationStateP& ss,
                                           PeridynamicsFlags* flags)
  : Uintah::Material(ps), d_materialModel(0), d_particle_creator(0)
{
  d_varLabel = scinew PeridynamicsLabel();

  // The standard set of initializations needed
  standardInitialization(ps, flags);
  
  d_materialModel->setSharedState(ss.get_rep());
  //d_failureModel->setSharedState(ss.get_rep());  // TODO: Not sure about this

  // Check to see which ParticleCreator object we need
  d_particle_creator = ParticleCreatorFactory::create(ps, this, flags);

  // Create and save a family computer instance
  d_family_computer = scinew FamilyComputer(flags, d_varLabel);
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

  d_failureModel = PeridynamicsFailureModelFactory::create(ps,flags);
  if(!d_failureModel){
    std::ostringstream desc;
    desc << "An error occured in the PeridynamicsFailureModelFactory that has \n" 
         << " slipped through the existing bullet proofing. " << std::endl;
    throw Uintah::ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  std::list<Uintah::GeometryObject::DataItem> geom_obj_data;
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
PeridynamicsMaterial::PeridynamicsMaterial() : d_materialModel(0), d_particle_creator(0), d_family_computer(0)
{
  d_varLabel = scinew PeridynamicsLabel();
  d_failureModel = 0;
}

PeridynamicsMaterial::~PeridynamicsMaterial()
{
  delete d_varLabel;
  delete d_materialModel;
  delete d_particle_creator;
  delete d_family_computer;
  delete d_failureModel;

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
  d_failureModel->outputProblemSpec(peridynamics_ps);
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

PeridynamicsFailureModel* 
PeridynamicsMaterial::getFailureModel() const
{
  return d_failureModel;
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

void 
PeridynamicsMaterial::createNeighborList(const Uintah::Patch* patch,
                                         Uintah::DataWarehouse* new_dw)
{
  d_family_computer->createNeighborList(this, patch, new_dw);
}

FamilyComputer* 
PeridynamicsMaterial::getFamilyComputer()
{
  return  d_family_computer;
}

