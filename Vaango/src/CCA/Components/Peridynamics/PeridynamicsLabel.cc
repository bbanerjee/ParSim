#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Grid/Variables/ParticleVariable.h>
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
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pDispGradLabel = Uintah::VarLabel::create("p.displacementGradient",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pDefGradLabel = Uintah::VarLabel::create("p.deformationGradient",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pStressLabel = Uintah::VarLabel::create("p.stress",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription() );
  pVolumeLabel = Uintah::VarLabel::create("p.volume",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pVolumeDeformedLabel = Uintah::VarLabel::create("p.volumedeformed",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pMassLabel = Uintah::VarLabel::create("p.mass",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pVelocityLabel = Uintah::VarLabel::create("p.velocity", 
			Uintah::ParticleVariable<Vector>::getTypeDescription() );
  pAccelerationLabel = Uintah::VarLabel::create("p.acceleration",
				   Uintah::ParticleVariable<Vector>::getTypeDescription()); 
  pExternalForceLabel = Uintah::VarLabel::create("p.externalforce",
			Uintah::ParticleVariable<Vector>::getTypeDescription() );
  pXLabel = Uintah::VarLabel::create("p.x",
			     Uintah::ParticleVariable<Point>::getTypeDescription(),
			     IntVector(0,0,0), Uintah::VarLabel::PositionVariable);
  pParticleIDLabel = Uintah::VarLabel::create("p.particleID",
			Uintah::ParticleVariable<long64>::getTypeDescription() );
  pVelGradLabel_preReloc = Uintah::VarLabel::create("p.velocityGradient+",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pDispGradLabel_preReloc = Uintah::VarLabel::create("p.displacementGradient+",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pDefGradLabel_preReloc = Uintah::VarLabel::create("p.deformationGradient+",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription());
  pStressLabel_preReloc = Uintah::VarLabel::create("p.stress+",
			Uintah::ParticleVariable<Matrix3>::getTypeDescription() );
  pVolumeLabel_preReloc = Uintah::VarLabel::create("p.volume+",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pMassLabel_preReloc = Uintah::VarLabel::create("p.mass+",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pVelocityLabel_preReloc = Uintah::VarLabel::create("p.velocity+", 
			Uintah::ParticleVariable<Vector>::getTypeDescription() );
  pAccelerationLabel_preReloc = Uintah::VarLabel::create("p.acceleration+",
				   Uintah::ParticleVariable<Vector>::getTypeDescription()); 
  pExternalForceLabel_preReloc = Uintah::VarLabel::create("p.externalforce+",
			Uintah::ParticleVariable<Vector>::getTypeDescription() );
  pXLabel_preReloc = Uintah::VarLabel::create( "p.x+",
			Uintah::ParticleVariable<Point>::getTypeDescription(),
			IntVector(0,0,0),
			Uintah::VarLabel::PositionVariable);
  pParticleIDLabel_preReloc = Uintah::VarLabel::create("p.particleID+",
			Uintah::ParticleVariable<long64>::getTypeDescription() );

  // Node Centered Variables
  gAccelerationLabel = Uintah::VarLabel::create("g.acceleration",
			Uintah::NCVariable<Vector>::getTypeDescription() );
  gMassLabel = Uintah::VarLabel::create("g.mass",
			Uintah::NCVariable<double>::getTypeDescription() );
  gVelocityLabel = Uintah::VarLabel::create("g.velocity",
			Uintah::NCVariable<Vector>::getTypeDescription() );
  gExternalForceLabel = Uintah::VarLabel::create( "g.externalforce",
			Uintah::NCVariable<Vector>::getTypeDescription() );
  gInternalForceLabel = Uintah::VarLabel::create( "g.internalforce",
			Uintah::NCVariable<Vector>::getTypeDescription() );
  gContactLabel       = Uintah::VarLabel::create( "g.contact",
			Uintah::NCVariable<int>::getTypeDescription() );
  gStressLabel   = Uintah::VarLabel::create( "g.stress",
                   Uintah::NCVariable<Matrix3>::getTypeDescription() );
  gVolumeLabel     = Uintah::VarLabel::create("g.volume",
			Uintah::NCVariable<double>::getTypeDescription());

  // Reduction variables
  partCountLabel = Uintah::VarLabel::create("particleCount",
				   sumlong_vartype::getTypeDescription());
  delTLabel = Uintah::VarLabel::create( "delT", delt_vartype::getTypeDescription() );
  StrainEnergyLabel = Uintah::VarLabel::create( "StrainEnergy",
			sum_vartype::getTypeDescription() );
  AccStrainEnergyLabel = Uintah::VarLabel::create( "AccStrainEnergy",
			max_vartype::getTypeDescription() );
  KineticEnergyLabel = Uintah::VarLabel::create( "KineticEnergy",
			sum_vartype::getTypeDescription() );
  ThermalEnergyLabel = Uintah::VarLabel::create( "ThermalEnergy",
			sum_vartype::getTypeDescription() );
  TotalMassLabel = Uintah::VarLabel::create( "TotalMass",
				 sum_vartype::getTypeDescription() );
  CenterOfMassPositionLabel = Uintah::VarLabel::create( "CenterOfMassPosition",
				 sumvec_vartype::getTypeDescription() );
  TotalMomentumLabel = Uintah::VarLabel::create( "TotalMomentum",
				 sumvec_vartype::getTypeDescription() );
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
  Uintah::VarLabel::destroy(pExternalForceLabel);
  Uintah::VarLabel::destroy(pExternalForceLabel_preReloc);
  Uintah::VarLabel::destroy(pXLabel);
  Uintah::VarLabel::destroy(pXLabel_preReloc);
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

  Uintah::VarLabel::destroy(AccStrainEnergyLabel);
  Uintah::VarLabel::destroy(StrainEnergyLabel);
  Uintah::VarLabel::destroy(KineticEnergyLabel);
  Uintah::VarLabel::destroy(TotalMassLabel);
  Uintah::VarLabel::destroy(CenterOfMassPositionLabel);
}
