#include <Commands/DepositIntoContainerResumeCommand.h>

using namespace dem;
void
DepositIntoContainerResumeCommand::execute(Assembly* assembly)
{
  assembly->deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());

  if (assembly->getMPIRank() == 0) {
    const Box& allContainer = assembly->getAllContainer();
    assembly->setContainer(Box(
      allContainer.getMinCorner().getX(), allContainer.getMinCorner().getY(),
      allContainer.getMinCorner().getZ(), allContainer.getMaxCorner().getX(),
      allContainer.getMaxCorner().getY(),
      dem::Parameter::getSingleton().parameter["trimHeight"]));
    assembly->buildBoundary(6, "trim_boundary_ini");
    char cstr[50];
    std::size_t endSnap = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["endSnap"]);
    assembly->trim(false, Assembly::combineString(cstr, "deposit_particle_", endSnap, 3),
         "trim_particle_ini");
  }
}
