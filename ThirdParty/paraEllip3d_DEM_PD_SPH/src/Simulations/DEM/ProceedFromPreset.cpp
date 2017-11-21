#include <Simulations/DEM/ProceedFromPreset.h>
#include <Boundary/BoundaryReader.h>
#include <InputOutput/InputParameter.h>

#include <iostream>
#include <string>

using namespace dem;
void
ProceedFromPreset::execute(DiscreteElements* dem)
{
  std::string boundaryFilename = util::getFilename("boundaryFilename");
  std::string particleFilename = util::getFilename("particleFilename");
  //std::cout << "Called ProceedFromPresetCommand with: "
  //          << " BoundaryFile = " << boundaryFilename << " and "
  //          << " ParticleFile = " << particleFilename << std::endl;

  dem->allowPatchDomainResize(Boundary::BoundaryID::ZPLUS);
  dem->deposit(boundaryFilename, particleFilename);
}
