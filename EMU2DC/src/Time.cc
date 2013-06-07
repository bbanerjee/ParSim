#include <Time.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Emu2DC;

Time::Time()
  : d_max_time(1.0), d_delT(1.0e-2), d_max_iter(100), d_cur_time(0.0)
{
}

Time::Time(const Uintah::ProblemSpecP& ps)
  : d_cur_time(0.0)
{
  initialize(ps);
}

Time::~Time()
{
}

void
Time::initialize(const Uintah::ProblemSpecP& ps)
{
  // get simulation time information
  Uintah::ProblemSpecP time_ps = ps->findBlock("Time");
  if (!time_ps) {
    throw Exception("**ERROR** <Time> tag not found", __FILE__, __LINE__); 
  }

  time_ps->require("max_time", d_max_time);
  time_ps->require("max_iterations", d_max_iter);
  time_ps->require("delt", d_delT);
  
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Time& time)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Max T = " << time.d_max_time << " del T = " << time.d_delT
        << " Max iter = " << time.d_max_iter << std::endl;
    return out;
  }

}
