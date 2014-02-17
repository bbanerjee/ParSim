#ifndef __VAANGO_PERIDYNAMICS_DOMAIN_BOUND_COND_H__
#define __VAANGO_PERIDYNAMICS_DOMAIN_BOUND_COND_H__

#include <Core/Geometry/Vector.h>
#include <Core/Grid/Variables/NCVariable.h>

namespace Vaango {

  class PeridynamicsDomainBoundCond {

  public:
    
    PeridynamicsDomainBoundCond();
    ~PeridynamicsDomainBoundCond();

    void setBoundaryCondition(const Uintah::Patch* patch,int dwi, const std::string& type,
                              Uintah::NCVariable<SCIRun::Vector>& variable,
                              std::string interp_type="linear");

    void setBoundaryCondition(const Uintah::Patch* patch,int dwi, const std::string& type,
                              Uintah::NCVariable<double>& variable,
                              std::string interp_type="linear");

  };

}

#endif
