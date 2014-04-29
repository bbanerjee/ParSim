#include <CCA/Components/Peridynamics/FamilyComputer/Bond.h>

#include <iostream>
#define _USE_MATH_DEFINE
#include <cmath>

using namespace Vaango;

Bond::Bond()
  :d_start(0),d_end(0),d_broken(false)
{
}

Bond::Bond(const Uintah::ParticleID& start, const Uintah::ParticleID& end)
  :d_start(start),d_end(end),d_broken(false)
{
}

Bond::~Bond()
{
}

bool 
Bond::operator==(const Bond& bond) const
{
  return ((d_start == bond.d_start && d_end == bond.d_end) ||
          (d_start == bond.d_end && d_end == bond.d_start));
}

namespace Vaango {

  std::ostream& operator<<(std::ostream& out, const Bond& bond)
  {
    out.setf(std::ios::floatfield);
    out.precision(3);
    out << "Bond: [" << bond.d_start << " - " << bond.d_end 
        << "], broken = " << std::boolalpha << bond.d_broken 
        << std::endl;
    return out;
  }
}
