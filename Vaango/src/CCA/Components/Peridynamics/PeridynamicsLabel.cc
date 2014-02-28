#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>

using namespace Vaango;

PeridynamicsLabel::PeridynamicsLabel()
{
  // Particle Variables
  //non PermanentParticleState
  pPressureLabel  = Uintah::VarLabel::create("p.pressure",
			Uintah::ParticleVariable<double>::getTypeDescription() );

  //PermanentParticleState
  pVelGradLabel = Uintah::VarLabel::create("p.velocityGradient",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pDispGradLabel = Uintah::VarLabel::create("p.displacementGradient",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pDefGradLabel = Uintah::VarLabel::create("p.deformationGradient",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pStressLabel = Uintah::VarLabel::create("p.stress",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pVolumeLabel = Uintah::VarLabel::create("p.volume",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pVolumeDeformedLabel = Uintah::VarLabel::create("p.volumedeformed",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pMassLabel = Uintah::VarLabel::create("p.mass",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pVelocityLabel = Uintah::VarLabel::create("p.velocity", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pDispLabel = Uintah::VarLabel::create("p.displacement", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pVelocityStarLabel = Uintah::VarLabel::create("p.velocitystar", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pAccelerationLabel = Uintah::VarLabel::create("p.acceleration",
				   Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription()); 
  pInternalForceLabel = Uintah::VarLabel::create("p.internalforce",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pExternalForceLabel = Uintah::VarLabel::create("p.externalforce",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pPositionLabel = Uintah::VarLabel::create("p.position",
			     Uintah::ParticleVariable<SCIRun::Point>::getTypeDescription(),
			     SCIRun::IntVector(0,0,0), Uintah::VarLabel::PositionVariable);
  pSizeLabel = Uintah::VarLabel::create("p.size",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pParticleIDLabel = Uintah::VarLabel::create("p.particleID",
			Uintah::ParticleVariable<Uintah::long64>::getTypeDescription() );
  pVelGradLabel_preReloc = Uintah::VarLabel::create("p.velocityGradient+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pDispGradLabel_preReloc = Uintah::VarLabel::create("p.displacementGradient+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pDefGradLabel_preReloc = Uintah::VarLabel::create("p.deformationGradient+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pStressLabel_preReloc = Uintah::VarLabel::create("p.stress+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pVolumeLabel_preReloc = Uintah::VarLabel::create("p.volume+",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pMassLabel_preReloc = Uintah::VarLabel::create("p.mass+",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pVelocityLabel_preReloc = Uintah::VarLabel::create("p.velocity+", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pDispLabel_preReloc = Uintah::VarLabel::create("p.displacement+", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pAccelerationLabel_preReloc = Uintah::VarLabel::create("p.acceleration+",
				   Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription()); 
  pInternalForceLabel_preReloc = Uintah::VarLabel::create("p.internalforce+",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pExternalForceLabel_preReloc = Uintah::VarLabel::create("p.externalforce+",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pPositionLabel_preReloc = Uintah::VarLabel::create( "p.position+",
			Uintah::ParticleVariable<SCIRun::Point>::getTypeDescription(),
			SCIRun::IntVector(0,0,0),
			Uintah::VarLabel::PositionVariable);
  pParticleIDLabel_preReloc = Uintah::VarLabel::create("p.particleID+",
			Uintah::ParticleVariable<Uintah::long64>::getTypeDescription() );
  pSizeLabel_preReloc = Uintah::VarLabel::create("p.size+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );

  // Node Centered Variables
  gAccelerationLabel = Uintah::VarLabel::create("g.acceleration",
			Uintah::NCVariable<SCIRun::Vector>::getTypeDescription() );
  gMassLabel = Uintah::VarLabel::create("g.mass",
			Uintah::NCVariable<double>::getTypeDescription() );
  gVelocityLabel = Uintah::VarLabel::create("g.velocity",
			Uintah::NCVariable<SCIRun::Vector>::getTypeDescription() );
  gVelocityStarLabel = Uintah::VarLabel::create("g.velocitystar",
			Uintah::NCVariable<SCIRun::Vector>::getTypeDescription() );
  gExternalForceLabel = Uintah::VarLabel::create( "g.externalforce",
			Uintah::NCVariable<SCIRun::Vector>::getTypeDescription() );
  gInternalForceLabel = Uintah::VarLabel::create( "g.internalforce",
			Uintah::NCVariable<SCIRun::Vector>::getTypeDescription() );
  gContactLabel       = Uintah::VarLabel::create( "g.contact",
			Uintah::NCVariable<int>::getTypeDescription() );
  gStressLabel   = Uintah::VarLabel::create( "g.stress",
                   Uintah::NCVariable<Uintah::Matrix3>::getTypeDescription() );
  gVolumeLabel     = Uintah::VarLabel::create("g.volume",
			Uintah::NCVariable<double>::getTypeDescription());

  pCellNAPIDLabel = Uintah::VarLabel::create("pc.NAPID",
			Uintah::CCVariable<short int>::getTypeDescription());

  // Reduction variables
  partCountLabel = Uintah::VarLabel::create("particleCount",
				   Uintah::sumlong_vartype::getTypeDescription());
  delTLabel = Uintah::VarLabel::create( "delT", Uintah::delt_vartype::getTypeDescription() );
  StrainEnergyLabel = Uintah::VarLabel::create( "StrainEnergy",
			Uintah::sum_vartype::getTypeDescription() );
  AccStrainEnergyLabel = Uintah::VarLabel::create( "AccStrainEnergy",
			Uintah::max_vartype::getTypeDescription() );
  KineticEnergyLabel = Uintah::VarLabel::create( "KineticEnergy",
			Uintah::sum_vartype::getTypeDescription() );
  CenterOfMassPositionLabel = Uintah::VarLabel::create( "CenterOfMassPosition",
				 Uintah::sumvec_vartype::getTypeDescription() );
  TotalMomentumLabel = Uintah::VarLabel::create( "TotalMomentum",
				 Uintah::sumvec_vartype::getTypeDescription() );
} 

PeridynamicsLabel::~PeridynamicsLabel()
{
  //non PermanentParticleState
  Uintah::VarLabel::destroy(pVolumeDeformedLabel);
  Uintah::VarLabel::destroy(pPressureLabel);

  //PermanentParticleState
  Uintah::VarLabel::destroy(pDefGradLabel);
  Uintah::VarLabel::destroy(pVelGradLabel);
  Uintah::VarLabel::destroy(pDispGradLabel);
  Uintah::VarLabel::destroy(pDefGradLabel_preReloc);
  Uintah::VarLabel::destroy(pVelGradLabel_preReloc);
  Uintah::VarLabel::destroy(pDispGradLabel_preReloc);
  Uintah::VarLabel::destroy(pStressLabel);
  Uintah::VarLabel::destroy(pStressLabel_preReloc);
  Uintah::VarLabel::destroy(pVolumeLabel);
  Uintah::VarLabel::destroy(pVolumeLabel_preReloc);
  Uintah::VarLabel::destroy(pMassLabel);
  Uintah::VarLabel::destroy(pMassLabel_preReloc);
  Uintah::VarLabel::destroy(pVelocityLabel);
  Uintah::VarLabel::destroy(pVelocityLabel_preReloc);
  Uintah::VarLabel::destroy(pDispLabel);
  Uintah::VarLabel::destroy(pDispLabel_preReloc);
  Uintah::VarLabel::destroy(pAccelerationLabel);
  Uintah::VarLabel::destroy(pAccelerationLabel_preReloc);
  Uintah::VarLabel::destroy(pInternalForceLabel);
  Uintah::VarLabel::destroy(pInternalForceLabel_preReloc);
  Uintah::VarLabel::destroy(pExternalForceLabel);
  Uintah::VarLabel::destroy(pExternalForceLabel_preReloc);
  Uintah::VarLabel::destroy(pPositionLabel);
  Uintah::VarLabel::destroy(pPositionLabel_preReloc);
  Uintah::VarLabel::destroy(pSizeLabel);
  Uintah::VarLabel::destroy(pSizeLabel_preReloc);
  Uintah::VarLabel::destroy(pParticleIDLabel);
  Uintah::VarLabel::destroy(pParticleIDLabel_preReloc);

  Uintah::VarLabel::destroy(gAccelerationLabel);
  Uintah::VarLabel::destroy(gMassLabel);
  Uintah::VarLabel::destroy(gVelocityLabel);
  Uintah::VarLabel::destroy(gExternalForceLabel);
  Uintah::VarLabel::destroy(gInternalForceLabel);
  Uintah::VarLabel::destroy(gContactLabel);
  Uintah::VarLabel::destroy(gStressLabel);
  Uintah::VarLabel::destroy(gVolumeLabel);

  Uintah::VarLabel::destroy(delTLabel);

  Uintah::VarLabel::destroy(pCellNAPIDLabel);

  Uintah::VarLabel::destroy(partCountLabel);
  Uintah::VarLabel::destroy(AccStrainEnergyLabel);
  Uintah::VarLabel::destroy(StrainEnergyLabel);
  Uintah::VarLabel::destroy(KineticEnergyLabel);
  Uintah::VarLabel::destroy(CenterOfMassPositionLabel);
}