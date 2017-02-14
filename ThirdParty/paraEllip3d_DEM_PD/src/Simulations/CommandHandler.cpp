#include <Simulations/CavityExpansion.h>
#include <Simulations/CavityExpansionResume.h>
#include <Simulations/Command.h>
#include <Simulations/CommandHandler.h>
#include <Simulations/CoupledFluidFlow.h>
#include <Simulations/DepositIntoContainer.h>
#include <Simulations/DepositIntoContainerResume.h>
#include <Simulations/IsotropicLoading.h>
#include <Simulations/OedometerLoading.h>
#include <Simulations/PeridynamicsPullOut.h>
#include <Simulations/PeridynamicsRigidInclusion.h>
#include <Simulations/PlaneStrainLoading.h>
#include <Simulations/ProceedFromPreset.h>
#include <Simulations/TriaxialLoading.h>
#include <Simulations/TrimParticles.h>
#include <Simulations/TrueTriaxialLoading.h>
#include <Simulations/TuneMassPercentage.h>

using namespace dem;

CommandP
CommandHandler::handleCommand(int simuType)
{

  switch (simuType) {
    case 001: // proceed from preset state
      return std::make_unique<ProceedFromPreset>();
      break;
    case 002: // tune mass-percentage from number-percentage on size
              // distribution curve by trial and error
      return std::make_unique<TuneMassPercentage>();
      break;
    case 003: // trim particles
      return std::make_unique<TrimParticles>();
      break;
    case 101: // deposit spatially scattered particles into a rigid container
      return std::make_unique<DepositIntoContainer>();
      break;
    case 102: // resume deposition using specified data file of particles and
              // boundaries
      return std::make_unique<DepositIntoContainerResume>();
      break;
    case 201: // isotropic type 1 - create an initial state with low confining
              // pressure
      return std::make_unique<IsotropicLoading>();
      break;
    case 202: // isotropic type 2 - increase confining pressure from sigmaInit
              // to sigmaEnd
      return std::make_unique<IsotropicLoading>();
      break;
    case 203: // isotropic type 3 - conduct loading-unloading-reloading path
      return std::make_unique<IsotropicLoading>();
      break;
    case 301: // odometer type 1 - increase loading pressure
      return std::make_unique<OedometerLoading>();
      break;
    case 302: // odometer type 2 - loading-unloading-reloading
      return std::make_unique<OedometerLoading>();
      break;
    case 401: // triaxial type 1 - constant confining pressure
      return std::make_unique<TriaxialLoading>();
      break;
    case 402: // triaxial type 2 - loading-unloading-reloading
      return std::make_unique<TriaxialLoading>();
      break;
    case 411: // plane strain type 1 - in x direction
      return std::make_unique<PlaneStrainLoading>();
      break;
    case 412: // plane strain type 2 - loading-unloading-reloading
      return std::make_unique<PlaneStrainLoading>();
      break;
    case 501: // true triaxial 1 - create confining stress state
      return std::make_unique<TrueTriaxialLoading>();
      break;
    case 502: // true triaxial 2 - increase stress in one direction
      return std::make_unique<TrueTriaxialLoading>();
      break;
    case 601: // expand particles inside a virtual cavity and see what occurs
      return std::make_unique<CavityExpansion>();
      break;
    case 602: // resume expanding particles inside a virtual cavity and see what
              // occurs
      return std::make_unique<CavityExpansionResume>();
      break;
    case 701: // couple with sonic fluid flow
      return std::make_unique<CoupledFluidFlow>();
      break;
    case 3001: // rigid Inclusion problem
      return std::make_unique<PeridynamicsRigidInclusion>();
      break;
    case 3002: // pull out DEM particle from center of peri-domain
      return std::make_unique<PeridynamicsPullOut>();
      break;
    default:
      return nullptr;
  }
  return nullptr;
}
