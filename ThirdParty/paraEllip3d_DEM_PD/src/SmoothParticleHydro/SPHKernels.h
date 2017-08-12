#ifndef SPH_KERNELS_H
#define SPH_KERNELS_H

#include <Core/Math/Vec.h>
#include <Core/Types/integertypes.h>
#include <Core/Types/realtypes.h>

#include <iostream>
#include <math.h>
#include <vector>

namespace sph {

class SPHKernels
{

  public:

  SPHKernels() = default;

  constexpr static const REAL qminQuinticSplineKernel = 0.7593;

  template <int dim>
  REAL minQuinticSplineKernel(const REAL& smoothLength) const;

  template <int dim>
  REAL quinticSplineKernel(const dem::Vec& pos_a,
                           const dem::Vec& pos_b,
                           const REAL& smoothLength) const;

  template <int dim>
  REAL quinticSplineKernel(const REAL& dist_ab,
                           const REAL& smoothLength) const;

  template <int dim>
  dem::Vec gradientQuinticSplineKernel(const dem::Vec& pos_a,
                                       const dem::Vec& pos_b,
                                       const REAL& smoothLength) const;

  template <int dim>
  dem::Vec gradientQuinticSplineKernel(const dem::Vec& vec_ab,
                                       const REAL& dist_ab,
                                       const REAL& smoothLength) const;

  private:

  template <int dim>
  REAL computeQuinticSplineKernelFactor(const REAL& smoothLength) const;

  REAL computeQuinticSplineKernel(const REAL& qval,
                                  const REAL& factor) const;

  REAL computeGradientQuinticSplineKernel(const REAL& qval,
                                          const REAL& factor) const;
}; 

} // end namespace sph

#endif
