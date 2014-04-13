#ifndef __MATITI_PROBLEMSPEC_UTIL_H__
#define __MATITI_PROBLEMSPEC_UTIL_H__

#include <Geometry/Polygon3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <string>

namespace Matiti_ProblemSpecUtil
{
  // A routine to read in the boundary of a two-dimensional region
  void readVector(Uintah::ProblemSpecP& ps, std::vector<double>& coeffVector);
  void readBoundary(Uintah::ProblemSpecP& ps, Matiti::Polygon3D& boundary);
  void parseVector(const std::string& stringValue, SCIRun::Vector& value);
  void checkForInputError(const std::string & stringValue);
}

#endif
