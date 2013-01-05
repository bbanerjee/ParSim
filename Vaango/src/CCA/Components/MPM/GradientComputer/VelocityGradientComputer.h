#ifndef __VELOCITY_GRADIENT_COMPUTER_H__
#define __VELOCITY_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>


namespace Uintah {

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class VelocityGradientComputer
    \brief Class for computing velocity gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class VelocityGradientComputer : public GradientComputer {

  public:
         
    VelocityGradientComputer(MPMFlags* MFlag);
    VelocityGradientComputer(const VelocityGradientComputer* gc);
    virtual ~VelocityGradientComputer();

    // Make a clone of the gradient computer
    VelocityGradientComputer* clone();

    // Actually compute velocity gradient
    void computeVelGrad(ParticleInterpolator* interpolator,
                        const double* oodx,
                        const short pgFld[],
                        const Point& px,
                        const Matrix3& pSize,
                        const Matrix3& pDefGrad_old,
                        constNCVariable<Vector>& gVelocity,
                        constNCVariable<Vector>& GVelocity,
                        Matrix3& velGrad_new);

  protected:
    void computeAxiSymVelocityGradient(Matrix3& velGrad,
                                       vector<IntVector>& ni,
                                       vector<Vector>& d_S,
                                       vector<double>& S,
                                       const double* oodx,
                                       constNCVariable<Vector>& gVelocity,
                                       const Point& px);

    void computeVelocityGradient(Matrix3& velGrad,
                                 vector<IntVector>& ni,
                                 vector<Vector>& d_S,
                                 const double* oodx, 
                                 const short pgFld[],
                                 constNCVariable<Vector>& gVelocity,
                                 constNCVariable<Vector>& GVelocity);
    
  };
} // End namespace Uintah
      


#endif  // __VELOCITY_GRADIENT_COMPUTER_H__

