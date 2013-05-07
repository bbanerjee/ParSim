#ifndef __EMU2DC_PROBLEMSPEC_UTIL_H__
#define __EMU2DC_PROBLEMSPEC_UTIL_H__

#include <Core/Geometry/Vector.h>
#include <string>

namespace Emu2DC_ProblemSpecUtil
{
  void parseVector(const std::string& stringValue, SCIRun::Vector& value);
  void checkForInputError(const std::string & stringValue);
}

#endif
