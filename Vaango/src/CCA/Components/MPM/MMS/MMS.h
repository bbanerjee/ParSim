/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2016 Parresai Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
/*
  MMS.cc -  Supports these manufactured solutions

  1) Axis Aligned MMS : Already was a part of Uintah. 
  Paper : An evaluation of explicit time integration schemes for use with the 
  generalized interpolation material point method ", Volume 227, pp.9628-9642 2008
  2) Generalized Vortex : Newly added
  Paper : Establishing Credibility of Particle Methods through Verification testing. 
  Particles 2011 II International Conference on Particle-based methods 
  Fundamentals and Applications.
  3) Expanding Ring : Newly added
  Paper : An evaluation of explicit time integration schemes for use with the 
  generalized interpolation material point method ", Volume 227, pp.9628-9642 2008
  4) Uniaxial-strain wave propagation:
  See ONR-MURI December 2015 monthly report PAR-10021867-1516.03, Appendix A.4.

  Member Functions :

  initializeParticleForMMS : Initializes the Particle data at t = 0 ; 
  Some MMS have intial velocity/displacement/stress. 
  For initial stress state, look at cnh_mms.cc 

  computeExternalForceForMMS : Computes the analytically determined body force for the 
  pre-determined deformation. 

  Author : Krishna Kamojjala
  Department of Mechanical Engineering
  University of Utah.
  Date   : 110824
        
*/

#ifndef __COMPONENTS_MPM_MMS_H__
#define __COMPONENTS_MPM_MMS_H__

#include <CCA/Components/MPM/SerialMPM.h>
#include <cmath>
#include <iostream>

namespace Uintah {

  class MMS {

  public :
        
    void initializeParticleForMMS(ParticleVariable<Point> &position,
                                  ParticleVariable<Vector> &pvelocity,
                                  ParticleVariable<Matrix3> &psize,
                                  ParticleVariable<Vector> &pdisp,
                                  ParticleVariable<double> &pmass,
                                  ParticleVariable<double> &pvolume ,
                                  Point p, 
                                  Vector dxcc, 
                                  Matrix3 size , 
                                  const Patch* patch,
                                  MPMFlags* flags,
                                  particleIndex i );

    void computeExternalForceForMMS(DataWarehouse* old_dw,
                                    DataWarehouse* new_dw, 
                                    double time, 
                                    ParticleSubset* pset, 
                                    MPMLabel* lb, 
                                    MPMFlags* flags , 
                                    ParticleVariable<Vector> &ExtForce);

  private :

    void initGeneralizedVortex(const MPMFlags* flags,
                               particleIndex pidx,
                               const Point& p,
                               const Vector& dxcc,
                               const Matrix3& size,
                               ParticleVariable<double>& pvolume,
                               ParticleVariable<double>& pmass,
                               ParticleVariable<Point>& position,
                               ParticleVariable<Vector>& pvelocity,
                               ParticleVariable<Vector>& pdisp,
                               ParticleVariable<Matrix3>& psize);

    void extForceGeneralizedVortex(const MPMFlags* flags,
                                   const MPMLabel* lb,
                                   const double& time,
                                   ParticleSubset* pset,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   ParticleVariable<Vector>& pExtForce);

    void initAxisAligned(const MPMFlags* flags,
                         particleIndex pidx,
                         const Point& p,
                         const Vector& dxcc,
                         const Matrix3& size,
                         ParticleVariable<double>& pvolume,
                         ParticleVariable<double>& pmass,
                         ParticleVariable<Point>& position,
                         ParticleVariable<Vector>& pvelocity,
                         ParticleVariable<Vector>& pdisp,
                         ParticleVariable<Matrix3>& psize);

    void extForceAxisAligned(const MPMFlags* flags,
                             const MPMLabel* lb,
                             const double& time,
                             ParticleSubset* pset,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             ParticleVariable<Vector>& pExtForce);

    void initExpandingRing(const MPMFlags* flags,
                           particleIndex pidx,
                           const Point& p,
                           const Vector& dxcc,
                           const Matrix3& size,
                           ParticleVariable<double>& pvolume,
                           ParticleVariable<double>& pmass,
                           ParticleVariable<Point>& position,
                           ParticleVariable<Vector>& pvelocity,
                           ParticleVariable<Vector>& pdisp,
                           ParticleVariable<Matrix3>& psize);

    void extForceExpandingRing(const MPMFlags* flags,
                               const MPMLabel* lb,
                               const double& time,
                               ParticleSubset* pset,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw,
                               ParticleVariable<Vector>& pExtForce);

    void initUniaxialStrain(const MPMFlags* flags,
                            particleIndex pidx,
                            const Point& p,
                            const Vector& dxcc,
                            const Matrix3& size,
                            ParticleVariable<double>& pvolume,
                            ParticleVariable<double>& pmass,
                            ParticleVariable<Point>& position,
                            ParticleVariable<Vector>& pvelocity,
                            ParticleVariable<Vector>& pdisp,
                            ParticleVariable<Matrix3>& psize);

    void bodyForceUniaxialStrainNonZeroInitStress(const MPMLabel* lb,
                                                  const double& time,
                                                  ParticleSubset* pset,
                                                  DataWarehouse* old_dw,
                                                  ParticleVariable<Vector>& pBodyForce);

    void extForceUniaxialStrainNonZeroInitStress(const MPMFlags* flags,
                                                 const MPMLabel* lb,
                                                 const double& time,
                                                 ParticleSubset* pset,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw,
                                                 ParticleVariable<Vector>& pExtForce);
  };

}// end namespace Uintah
#endif
