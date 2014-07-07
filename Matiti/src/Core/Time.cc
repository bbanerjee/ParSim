#include <Core/Time.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Matiti;

Time::Time()
  : d_max_time(1.0), d_delT(1.0e-2), d_max_iter(100), d_factor(0.5), d_cur_time(0.0)
{
}

Time::Time(double maxTime, double delT, int maxIter, double factor,
           double curTime)
  : d_max_time(maxTime), d_delT(delT), d_max_iter(maxIter),
    d_factor(factor), d_cur_time(curTime)
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
Time::clone(const Time& time)
{
  d_max_time = time.maxTime();
  d_delT = time.delT(); 
  d_max_iter = time.maxIter();
  d_factor = time.timeStepFactor(); 
  d_cur_time = time.currentTime();
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
  d_factor = 0.5;
  time_ps->get("time_step_reduction_factor", d_factor);
  
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Time& time)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "t = " << time.d_cur_time << " max t = " << time.d_max_time << " del t = " << time.d_delT
        << " max iter = " << time.d_max_iter << " factor = " << time.d_factor << std::endl;
    return out;
  }

}
