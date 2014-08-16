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
  pVelocityStarLabel = Uintah::VarLabel::create("p.velocitystar", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pAccelerationLabel = Uintah::VarLabel::create("p.acceleration",
				   Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription()); 
  pInternalForceLabel = Uintah::VarLabel::create("p.internalforce",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pExternalForceLabel = Uintah::VarLabel::create("p.externalforce",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pSizeLabel = Uintah::VarLabel::create("p.size",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pParticleIDLabel = Uintah::VarLabel::create("p.particleID",
			Uintah::ParticleVariable<Uintah::long64>::getTypeDescription() );

  // Labels for moving around within patches
  pStressLabel_preReloc = Uintah::VarLabel::create("p.stress+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pVolumeLabel_preReloc = Uintah::VarLabel::create("p.volume+",
			Uintah::ParticleVariable<double>::getTypeDescription());
  pMassLabel_preReloc = Uintah::VarLabel::create("p.mass+",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pVelocityLabel_preReloc = Uintah::VarLabel::create("p.velocity+", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pAccelerationLabel_preReloc = Uintah::VarLabel::create("p.acceleration+",
				   Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription()); 
  pInternalForceLabel_preReloc = Uintah::VarLabel::create("p.internalforce+",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pExternalForceLabel_preReloc = Uintah::VarLabel::create("p.externalforce+",
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
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

  // MPM Physical BC labels (permanent particle state)
  surfaceParticlesPerLoadCurveLabel = Uintah::VarLabel::create("particlesPerCurve",
                            Uintah::sumlong_vartype::getTypeDescription());
  pLoadCurveIDLabel = Uintah::VarLabel::create("p.loadCurveID",
                            Uintah::ParticleVariable<int>::getTypeDescription());
  pLoadCurveIDLabel_preReloc = Uintah::VarLabel::create("p.loadCurveID+",
                            Uintah::ParticleVariable<int>::getTypeDescription());

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

  // Peridynamics labels
  pPositionLabel = Uintah::VarLabel::create("p.x",
			     Uintah::ParticleVariable<SCIRun::Point>::getTypeDescription(),
			     SCIRun::IntVector(0,0,0), Uintah::VarLabel::PositionVariable);
  pHorizonLabel = Uintah::VarLabel::create("p.horizon",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pDamageLabel = Uintah::VarLabel::create("p.damage",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pDisplacementLabel = Uintah::VarLabel::create("p.displacement", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pDefGradLabel = Uintah::VarLabel::create("p.deformationGradient",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pShapeTensorInvLabel = Uintah::VarLabel::create("p.shapeTensorInverse",
                        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pPK1StressLabel = Uintah::VarLabel::create("p.PK1stress",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pNeighborListLabel = Uintah::VarLabel::create("p.neighborlist",
			Uintah::ParticleVariable<Uintah::NeighborList>::getTypeDescription() );
  pNeighborConnLabel = Uintah::VarLabel::create("p.neighborconn",
	                Uintah::ParticleVariable<Uintah::NeighborConnectivity>::getTypeDescription() );
  pNeighborCountLabel =  Uintah::VarLabel::create("p.neighborcount",
			Uintah::ParticleVariable<int>::getTypeDescription() );
  pNeighborBondEnergyLabel =  Uintah::VarLabel::create("p.bondEnergy",
			Uintah::ParticleVariable<Uintah::NeighborBondEnergy>::getTypeDescription() );
  pNeighborBondForceLabel =  Uintah::VarLabel::create("p.bondForce",
			Uintah::ParticleVariable<Uintah::NeighborBondInternalForce>::getTypeDescription() );

  pPositionStarLabel = Uintah::VarLabel::create( "p.positionstar",
			Uintah::ParticleVariable<SCIRun::Point>::getTypeDescription() );
  pDisplacementStarLabel = Uintah::VarLabel::create("p.displacementstar", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );

  pPositionLabel_preReloc = Uintah::VarLabel::create( "p.x+",
			Uintah::ParticleVariable<SCIRun::Point>::getTypeDescription(),
			SCIRun::IntVector(0,0,0),
			Uintah::VarLabel::PositionVariable);
  pHorizonLabel_preReloc = Uintah::VarLabel::create("p.horizon+",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pDamageLabel_preReloc = Uintah::VarLabel::create("p.damage+",
			Uintah::ParticleVariable<double>::getTypeDescription() );
  pDisplacementLabel_preReloc = Uintah::VarLabel::create("p.displacement+", 
			Uintah::ParticleVariable<SCIRun::Vector>::getTypeDescription() );
  pDefGradLabel_preReloc = Uintah::VarLabel::create("p.deformationGradient+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pShapeTensorInvLabel_preReloc = Uintah::VarLabel::create("p.shapeTensorInverse+",
                        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pPK1StressLabel_preReloc = Uintah::VarLabel::create("p.PK1stress+",
			Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription() );
  pNeighborListLabel_preReloc = Uintah::VarLabel::create("p.neighborlist+",
			Uintah::ParticleVariable<Uintah::NeighborList>::getTypeDescription() );
  pNeighborConnLabel_preReloc = Uintah::VarLabel::create("p.neighborconn+",
	                Uintah::ParticleVariable<Uintah::NeighborConnectivity>::getTypeDescription() );
  pNeighborCountLabel_preReloc =  Uintah::VarLabel::create("p.neighborcount+",
			Uintah::ParticleVariable<int>::getTypeDescription() );
  pNeighborBondEnergyLabel_preReloc =  Uintah::VarLabel::create("p.bondEnergy+",
			Uintah::ParticleVariable<Uintah::NeighborBondEnergy>::getTypeDescription() );
  pNeighborBondForceLabel_preReloc =  Uintah::VarLabel::create("p.bondForce+",
			Uintah::ParticleVariable<Uintah::NeighborBondInternalForce>::getTypeDescription() );


} 

PeridynamicsLabel::~PeridynamicsLabel()
{
  //non PermanentParticleState
  Uintah::VarLabel::destroy(pVolumeDeformedLabel);
  Uintah::VarLabel::destroy(pPressureLabel);

  //PermanentParticleState
  Uintah::VarLabel::destroy(pDefGradLabel);
  Uintah::VarLabel::destroy(pDefGradLabel_preReloc);
  Uintah::VarLabel::destroy(pStressLabel);
  Uintah::VarLabel::destroy(pStressLabel_preReloc);
  Uintah::VarLabel::destroy(pVolumeLabel);
  Uintah::VarLabel::destroy(pVolumeLabel_preReloc);
  Uintah::VarLabel::destroy(pMassLabel);
  Uintah::VarLabel::destroy(pMassLabel_preReloc);
  Uintah::VarLabel::destroy(pVelocityLabel);
  Uintah::VarLabel::destroy(pVelocityLabel_preReloc);
  Uintah::VarLabel::destroy(pAccelerationLabel);
  Uintah::VarLabel::destroy(pAccelerationLabel_preReloc);
  Uintah::VarLabel::destroy(pInternalForceLabel);
  Uintah::VarLabel::destroy(pInternalForceLabel_preReloc);
  Uintah::VarLabel::destroy(pExternalForceLabel);
  Uintah::VarLabel::destroy(pExternalForceLabel_preReloc);
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

  //--------------------------------------
  // Peridynamics labels
  //--------------------------------------
  Uintah::VarLabel::destroy(pPositionLabel);
  Uintah::VarLabel::destroy(pHorizonLabel);
  Uintah::VarLabel::destroy(pDamageLabel);
  Uintah::VarLabel::destroy(pDisplacementLabel);
  Uintah::VarLabel::destroy(pPK1StressLabel);
  Uintah::VarLabel::destroy(pNeighborListLabel);
  Uintah::VarLabel::destroy(pNeighborConnLabel);
  Uintah::VarLabel::destroy(pNeighborCountLabel);
  Uintah::VarLabel::destroy(pNeighborBondEnergyLabel);

  Uintah::VarLabel::destroy(pPositionStarLabel);
  Uintah::VarLabel::destroy(pDisplacementStarLabel);

  Uintah::VarLabel::destroy(pPositionLabel_preReloc);
  Uintah::VarLabel::destroy(pHorizonLabel_preReloc);
  Uintah::VarLabel::destroy(pDamageLabel_preReloc);
  Uintah::VarLabel::destroy(pDisplacementLabel_preReloc);
  Uintah::VarLabel::destroy(pPK1StressLabel_preReloc);
  Uintah::VarLabel::destroy(pNeighborListLabel_preReloc);
  Uintah::VarLabel::destroy(pNeighborConnLabel_preReloc);
  Uintah::VarLabel::destroy(pNeighborCountLabel_preReloc);
  Uintah::VarLabel::destroy(pNeighborBondEnergyLabel_preReloc);
  //--------------------------------------

}
