#include <Commands/Command.h>
#include <Commands/CommandHandler.h>
#include <Commands/CavityExpansionCommand.h>
#include <Commands/CavityExpansionResumeCommand.h>
#include <Commands/CoupledFluidFlowCommand.h>
#include <Commands/DepositIntoContainerCommand.h>
#include <Commands/DepositIntoContainerResumeCommand.h>
#include <Commands/IsotropicLoadingCommand.h>
#include <Commands/OedometerLoadingCommand.h>
#include <Commands/PeridynamicsPullOutCommand.h>
#include <Commands/PeridynamicsRigidInclusionCommand.h>
#include <Commands/PlaneStrainLoadingCommand.h>
#include <Commands/ProceedFromPresetCommand.h>
#include <Commands/TriaxialLoadingCommand.h>
#include <Commands/TrimParticlesCommand.h>
#include <Commands/TrueTriaxialLoadingCommand.h>
#include <Commands/TuneMassPercentageCommand.h>


using namespace dem;

CommandP
CommandHandler::handleCommand(int simuType)
{

  switch (simuType) {
    case 001: // proceed from preset state
      return std::make_unique<ProceedFromPresetCommand>();
      break;
    case 002: // tune mass-percentage from number-percentage on size
              // distribution curve by trial and error
      return std::make_unique<TuneMassPercentageCommand>();
      break;
    case 003: // trim particles
      return std::make_unique<TrimParticlesCommand>();
      break;
    case 101: // deposit spatially scattered particles into a rigid container
      return std::make_unique<DepositIntoContainerCommand>();
      break;
    case 102: // resume deposition using specified data file of particles and
              // boundaries
      return std::make_unique<DepositIntoContainerResumeCommand>();
      break;
    case 201: // isotropic type 1 - create an initial state with low confining
              // pressure
      return std::make_unique<IsotropicLoadingCommand>();
      break;
    case 202: // isotropic type 2 - increase confining pressure from sigmaInit
              // to sigmaEnd
      return std::make_unique<IsotropicLoadingCommand>();
      break;
    case 203: // isotropic type 3 - conduct loading-unloading-reloading path
      return std::make_unique<IsotropicLoadingCommand>();
      break;
    case 301: // odometer type 1 - increase loading pressure
      return std::make_unique<OedometerLoadingCommand>();
      break;
    case 302: // odometer type 2 - loading-unloading-reloading
      return std::make_unique<OedometerLoadingCommand>();
      break;
    case 401: // triaxial type 1 - constant confining pressure
      return std::make_unique<TriaxialLoadingCommand>();
      break;
    case 402: // triaxial type 2 - loading-unloading-reloading
      return std::make_unique<TriaxialLoadingCommand>();
      break;
    case 411: // plane strain type 1 - in x direction
      return std::make_unique<PlaneStrainLoadingCommand>();
      break;
    case 412: // plane strain type 2 - loading-unloading-reloading
      return std::make_unique<PlaneStrainLoadingCommand>();
      break;
    case 501: // true triaxial 1 - create confining stress state
      return std::make_unique<TrueTriaxialLoadingCommand>();
      break;
    case 502: // true triaxial 2 - increase stress in one direction
      return std::make_unique<TrueTriaxialLoadingCommand>();
      break;
    case 601: // expand particles inside a virtual cavity and see what occurs
      return std::make_unique<CavityExpansionCommand>();
      break;
    case 602: // resume expanding particles inside a virtual cavity and see what
              // occurs
      return std::make_unique<CavityExpansionResumeCommand>();
      break;
    case 701: // couple with sonic fluid flow
      return std::make_unique<CoupledFluidFlowCommand>();
      break;
    case 3001: // rigid Inclusion problem
      return std::make_unique<PeridynamicsRigidInclusionCommand>();
      break;
    case 3002: // pull out DEM particle from center of peri-domain
      return std::make_unique<PeridynamicsPullOutCommand>();
      break;
    default:
      return nullptr;
  }
  return nullptr;
}