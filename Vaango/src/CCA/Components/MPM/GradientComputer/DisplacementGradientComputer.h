#ifndef __DISPLACEMENT_GRADIENT_COMPUTER_H__
#define __DISPLACEMENT_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>


namespace Uintah {

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class DisplacementGradientComputer
    \brief Class for computing displacement gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class DisplacementGradientComputer : public GradientComputer {

  public:
         
    DisplacementGradientComputer(MPMFlags* MFlag);
    DisplacementGradientComputer(const DisplacementGradientComputer* gc);
    virtual ~DisplacementGradientComputer();

    // Make a clone of the gradient computer
    DisplacementGradientComputer* clone();

    // Actually compute displacement gradient
    void computeDispGrad(ParticleInterpolator* interp,
                         const double* oodx,
                         const Point& px,
                         const Matrix3& psize,
                         const Matrix3& pDefGrad_old,
                         constNCVariable<Vector> gDisp,
                         Matrix3& dispGrad_new);

  };
} // End namespace Uintah
      


#endif  // __DISPLACEMENT_GRADIENT_COMPUTER_H__

