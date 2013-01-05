#ifndef __GRADIENT_COMPUTER_H__
#define __GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/LinearInterpolator.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <vector>


namespace Uintah {

  class MPMFlags;
  class ParticleVariableBase;

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class GradientComputer
    \brief Base class for computing gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class GradientComputer {
  public:
         
    GradientComputer(MPMFlags* MFlag);
    GradientComputer(const GradientComputer* gc);
    virtual ~GradientComputer();

    // Make a clone of the gradient computer
    virtual GradientComputer* clone() = 0;

  public:

    /*! Calculate gradient of vector field for 8 noded interpolation, B matrix
        for Kmat and B matrix for Kgeo */
    void computeGradAndBmats(Matrix3& grad,
                             vector<IntVector>& ni,
                             vector<Vector>& d_S,
                             const double* oodx, 
                             constNCVariable<Vector>& gVec,
                             const Array3<int>& l2g,
                             double B[6][24],
                             double Bnl[3][24],
                             int* dof);

    void computeBmats(vector<IntVector>& ni,
                      vector<Vector>& d_S,
                      const double* oodx, 
                      const Array3<int>& l2g,
                      double B[6][24],
                      double Bnl[3][24],
                      int* dof);

  protected:
    
    /*! Calculate gradient of a vector field for 8 noded interpolation */
    void computeGrad(Matrix3& grad,
                     vector<IntVector>& ni,
                     vector<Vector>& d_S,
                     const double* oodx, 
                     constNCVariable<Vector>& gVec);

    MPMFlags* flag;

  };
} // End namespace Uintah
      

#endif  // __GRADIENT_COMPUTER_H__

