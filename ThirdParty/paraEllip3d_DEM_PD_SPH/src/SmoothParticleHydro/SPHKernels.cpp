#include <SmoothParticleHydro/SPHKernels.h>
#include <Core/Const/Constants.h>

using namespace sph;

using dem::Vec;

template <int dim>
REAL 
SPHKernels::minQuinticSplineKernel(const REAL& smoothLength) const
{
  REAL qval = SPHKernels::qminQuinticSplineKernel;
  REAL factor = computeQuinticSplineKernelFactor<dim>(smoothLength);
  return computeQuinticSplineKernel(qval, factor);
}

template <int dim>
REAL 
SPHKernels::quinticSplineKernel(const dem::Vec& pos_a,
                                const dem::Vec& pos_b,
                                const REAL& smoothLength) const
{
  Vec vec_ab = pos_a - pos_b;
  REAL dist_ab = vec_ab.length();
  return quinticSplineKernel<dim>(dist_ab, smoothLength);
}

template <int dim>
REAL 
SPHKernels::quinticSplineKernel(const REAL& dist_ab,
                                const REAL& smoothLength) const
{
  REAL qval = dist_ab/smoothLength;
  REAL factor = computeQuinticSplineKernelFactor<dim>(smoothLength);
  return computeQuinticSplineKernel(qval, factor);
}

template <int dim>
dem::Vec 
SPHKernels::gradientQuinticSplineKernel(const dem::Vec& pos_a,
                                        const dem::Vec& pos_b,
                                        const REAL& smoothLength) const
{
  dem::Vec vec_ab = pos_a - pos_b;
  REAL dist_ab = vec_ab.length();
  return gradientQuinticSplineKernel<dim>(vec_ab, dist_ab, smoothLength);
}

template <int dim>
dem::Vec 
SPHKernels::gradientQuinticSplineKernel(const dem::Vec& vec_ab,
                                        const REAL& dist_ab,
                                        const REAL& smoothLength) const
{
  REAL qval = dist_ab/smoothLength;
  REAL factor = computeQuinticSplineKernelFactor<dim>(smoothLength);
  REAL gradientFactor = factor/smoothLength;
  REAL gradient = 0.0;
  if (dist_ab > 1.0e-12) {
    gradient = computeGradientQuinticSplineKernel(qval, gradientFactor);
  }
  return (vec_ab/dist_ab)*gradient;
}

template <>
REAL 
SPHKernels::computeQuinticSplineKernelFactor<1>(const REAL& smoothLength) const
{
  return 1.0/(120.0*dem::Pi*smoothLength); 
}

template <>
REAL 
SPHKernels::computeQuinticSplineKernelFactor<2>(const REAL& smoothLength) const
{
  return 7.0/(478.0*dem::Pi*smoothLength*smoothLength); 
}

template <>
REAL 
SPHKernels::computeQuinticSplineKernelFactor<3>(const REAL& smoothLength) const
{
  return 3.0/(359.0*dem::Pi*smoothLength*smoothLength*smoothLength); 
}

REAL 
SPHKernels::computeQuinticSplineKernel(const REAL& qval,
                                       const REAL& factor) const
{
  REAL kernel = 0.0;
  REAL fac1 = 1.0 - qval;
  REAL fac2 = 2.0 - qval;
  REAL fac3 = 3.0 - qval;
  if (qval > 3.0) {
    kernel = 0.0;
  } else if (qval > 2.0) {
    REAL term1 = fac3*fac3*fac3*fac3*fac3;
    kernel =  factor * term1;
  } else if (qval > 1.0) {
    REAL term1 = fac3*fac3*fac3*fac3*fac3;
    REAL term2 = fac2*fac2*fac2*fac2*fac2;
    kernel =  factor * (term1 - 6 * term2);
  } else {
    REAL term1 = fac3*fac3*fac3*fac3*fac3;
    REAL term2 = fac2*fac2*fac2*fac2*fac2;
    REAL term3 = fac1*fac1*fac1*fac1*fac1;
    kernel =  factor * (term1 - 6 * term2 + 15 * term3);
  } 
  return kernel;
}

REAL
SPHKernels::computeGradientQuinticSplineKernel(const REAL& qval,
                                               const REAL& factor) const
{
  REAL gradient = 0.0;
  REAL fac1 = 1.0 - qval;
  REAL fac2 = 2.0 - qval;
  REAL fac3 = 3.0 - qval;

  if (qval > 3.0) {
    gradient = 0.0;
  } else if (qval > 2.0) {
    REAL term1 = fac3*fac3*fac3*fac3;
    gradient = factor * (-5.0 * term1);
  } else if (qval > 1.0) {
    REAL term1 = fac3*fac3*fac3*fac3;
    REAL term2 = fac2*fac2*fac2*fac2;
    gradient =  factor * (-5.0 * term1 + 30.0 * term2);
  } else {
    REAL term1 = fac3*fac3*fac3*fac3;
    REAL term2 = fac2*fac2*fac2*fac2;
    REAL term3 = fac1*fac1*fac1*fac1;
    gradient =  factor * (-5.0 * term1 + 30.0 * term2 - 75.0 * term3);
  } 
  return gradient;
}

namespace sph {

template REAL 
SPHKernels::minQuinticSplineKernel<1>(const REAL& smoothLength) const;
template REAL 
SPHKernels::minQuinticSplineKernel<2>(const REAL& smoothLength) const;
template REAL 
SPHKernels::minQuinticSplineKernel<3>(const REAL& smoothLength) const;

template REAL 
SPHKernels::quinticSplineKernel<1>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;
template REAL 
SPHKernels::quinticSplineKernel<2>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;
template REAL 
SPHKernels::quinticSplineKernel<3>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;

template REAL 
SPHKernels::quinticSplineKernel<1>(const REAL& dist_ab,
  const REAL& smoothLength) const;
template REAL 
SPHKernels::quinticSplineKernel<2>(const REAL& dist_ab,
  const REAL& smoothLength) const;
template REAL 
SPHKernels::quinticSplineKernel<3>(const REAL& dist_ab,
  const REAL& smoothLength) const;

template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<1>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;
template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<2>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;
template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<3>(const dem::Vec& pos_a,
  const dem::Vec& pos_b, const REAL& smoothLength) const;

template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<1>(const dem::Vec& vec_ab,
  const REAL& dist_ab, const REAL& smoothLength) const;
template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<2>(const dem::Vec& vec_ab,
  const REAL& dist_ab, const REAL& smoothLength) const;
template dem::Vec 
SPHKernels::gradientQuinticSplineKernel<3>(const dem::Vec& vec_ab,
  const REAL& dist_ab, const REAL& smoothLength) const;
}