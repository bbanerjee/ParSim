#ifndef __EMU2DC_PROBLEMSPEC_UTIL_H__
#define __EMU2DC_PROBLEMSPEC_UTIL_H__

#include <Geometry/Polygon3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <string>

namespace Emu2DC_ProblemSpecUtil
{
  // A routine to read in the boundary of a two-dimensional region
  void readBoundary(Uintah::ProblemSpecP& ps, Emu2DC::Polygon3D& boundary);
  void parseVector(const std::string& stringValue, SCIRun::Vector& value);
  void checkForInputError(const std::string & stringValue);
}

#endif
