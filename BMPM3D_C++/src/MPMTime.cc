#include <MPMTime.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

using namespace BrMPM;

MPMTime::MPMTime()
  : d_init_time(0.0), d_final_time(1.0), d_cfl(0.4), Time()
{
}

MPMTime::MPMTime(const Uintah::ProblemSpecP& ps)
  : d_cur_time(0.0), Time(ps)
{
}

MPMTime::~MPMTime()
{
}

void
MPMTime::initialize(const Uintah::ProblemSpecP& ps)
{
  // get simulation time information
  Uintah::ProblemSpecP time_ps = ps->findBlock("Time");
  if (!time_ps) {
    throw Exception("**ERROR** <Time> tag not found", __FILE__, __LINE__); 
  }

  time_ps->require("init_time", d_init_time);
  time_ps->require("final_time", d_final_time);
  time_ps->require("max_iterations", d_max_iter);
  time_ps->require("CFL", d_cfl);
  d_factor = 0.5;
  time_ps->get("time_step_reduction_factor", d_factor);
//  d_delT = std::min(d_domain.xrange()/(d_domain.numCells())[0], d_domain.yrange()/(d_domain.numCells())[1])
//           *d_cfl/(std::sqrt(d_mat->density()*d_mat->youngModulus());
}

/*namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const TTime& time)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "t = " << time.d_cur_time << " final t = " << time.d_final_time << " del t = " << time.d_delT
        << " max iter = " << time.d_max_iter << " factor = " << time.d_factor << std::endl;
    return out;
  }   

}  */
