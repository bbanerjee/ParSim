/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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
#include <CCA/Components/MPM/MMS/MMS.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Uintah;


//=====================================================================================
// Initialize MMS models
//=====================================================================================
void
MMS::initializeParticleForMMS(ParticleVariable<Point> &position,
                              ParticleVariable<Vector> &pvelocity,
                              ParticleVariable<Matrix3> &pSize,
                              ParticleVariable<Vector> &pdisp,
                              ParticleVariable<double> &pmass,
                              ParticleVariable<double> &pvolume ,
                              Point p, 
                              Vector dxcc, 
                              Matrix3 size , 
                              const Patch* patch,
                              MPMFlags* flags,
                              particleIndex pidx)
{

  string mms_type = flags->d_mmsType;
  if (mms_type == "GeneralizedVortex") {

    initGeneralizedVortex(flags, pidx, p, dxcc, size, 
                          pvolume, pmass, position, pvelocity, pdisp, pSize);

  } else if (mms_type == "AxisAligned" || mms_type == "AxisAligned3L") {

    initAxisAligned(flags, pidx, p, dxcc, size, 
                    pvolume, pmass, position, pvelocity, pdisp, pSize);

  } else if (mms_type == "ExpandingRing") {

    initExpandingRing(flags, pidx, p, dxcc, size, 
                      pvolume, pmass, position, pvelocity, pdisp, pSize);

  } else if (mms_type == "UniaxialStrainHarmonic") {

    initUniaxialStrainHarmonic(flags, pidx, p, dxcc, size, 
                               pvolume, pmass, position, pvelocity, pdisp, pSize);

  } else if (mms_type == "UniaxialStrainHomogeneousLinear") {

    initUniaxialStrainHomogeneousLinear(flags, pidx, p, dxcc, size, 
                                        pvolume, pmass, position, pvelocity, pdisp, pSize);

  } else if (mms_type == "UniaxialStrainHomogeneousQuadratic") {

    initUniaxialStrainHomogeneousQuadratic(flags, pidx, p, dxcc, size, 
                                           pvolume, pmass, position, pvelocity, pdisp, pSize);

  }
}

//=====================================================================================
// Compute external force for MMS models
//=====================================================================================
void
MMS::computeExternalForceForMMS(DataWarehouse* old_dw,
                                DataWarehouse* new_dw,
                                double time,
                                ParticleSubset* pset, 
                                MPMLabel* lb, 
                                MPMFlags* flags,
                                ParticleVariable<Vector> &pExtForce)
{   
  string mms_type = flags->d_mmsType;
  if (mms_type == "GeneralizedVortex") {

    extForceGeneralizedVortex(flags, lb,
                              time, pset, old_dw, new_dw,
                              pExtForce);

  } else if (mms_type == "AxisAligned" || mms_type == "AxisAligned3L") {

    extForceAxisAligned(flags, lb,
                        time, pset, old_dw, new_dw,
                        pExtForce);

  } else if (mms_type == "ExpandingRing") {

    extForceExpandingRing(flags, lb,
                          time, pset, old_dw, new_dw,
                          pExtForce);

  }
}

//=====================================================================================
// Compute body force for MMS models
//=====================================================================================
void
MMS::computeBodyForceForMMS(DataWarehouse* old_dw,
                            DataWarehouse* new_dw,
                            double time,
                            ParticleSubset* pset, 
                            MPMLabel* lb, 
                            MPMFlags* flags,
                            ParticleVariable<Vector> &pBodyForce)
{   
  string mms_type = flags->d_mmsType;

  if (mms_type == "UniaxialStrainHarmonic") {

    bodyForceUniaxialStrainHarmonic(lb, time, pset, old_dw, pBodyForce);

  } else if (mms_type == "UniaxialStrainHomogeneousLinear") {

    bodyForceUniaxialStrainHomogeneousLinear(lb, time, pset, old_dw, pBodyForce);

  } else if (mms_type == "UniaxialStrainHomogeneousQuadratic") {

    bodyForceUniaxialStrainHomogeneousQuadratic(lb, time, pset, old_dw, pBodyForce);

  }

}

//=====================================================================================
// Generalized vortex
//=====================================================================================
void 
MMS::initGeneralizedVortex(const MPMFlags* flags,
                           particleIndex pidx,
                           const Point& p,
                           const Vector& dxcc,
                           const Matrix3& size,
                           ParticleVariable<double>& pvolume,
                           ParticleVariable<double>& pmass,
                           ParticleVariable<Point>& position,
                           ParticleVariable<Vector>& pvelocity,
                           ParticleVariable<Vector>& pdisp,
                           ParticleVariable<Matrix3>& pSize)
{
  double t = 0.0;
  double A = 1.0;
  double R = sqrt(p.x()*p.x() + p.y()*p.y());
  double alfap = M_PI*(1. - 32.*pow((R-1),2.) + 256.*pow((R-1),4.));
  double alpha = A*sin(M_PI*t)*(1. - 32.*pow((R-1),2.) + 256.*pow((R-1),4.));
  double u = alfap*(-sin(alpha)*p.x() - cos(alpha)*p.y());
  double v = alfap*(cos(alpha)*p.x() - sin(alpha)*p.y());
  double w = 0.0;
  double rho0 = 1000.0;
  position[pidx] = p;
  if (flags->d_axisymmetric) {
    // assume unit radian extent in the circumferential direction
    pvolume[pidx]  = p.x()*(size(0,0)*size(1,1)-size(0,1)*size(1,0))*dxcc.x()*dxcc.y();
  } else {
    // standard voxel volume
    pvolume[pidx]  = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  }
  pvelocity[pidx]    = Vector(u,v,w);
  pSize[pidx]        = size;
  pmass[pidx]        = rho0*pvolume[pidx];
  pdisp[pidx]        = Vector(0.,0.,0.);
}

void 
MMS::extForceGeneralizedVortex(const MPMFlags* flags,
                               const MPMLabel* lb,
                               const double& time,
                               ParticleSubset* pset,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw,
                               ParticleVariable<Vector>& pExtForce)
{
  // std::cout << "Entered the GV loop " << std::endl;
  double mu = 384.615;
  double E = 1000.0;
  double rho0=1000.0;
  double c = sqrt(E/rho0);
  double A = 1.0;

  constParticleVariable<Point>  px;
  constParticleVariable<Vector> pdisp;
  constParticleVariable<double> pmass;
  old_dw->get(px,             lb->pXLabel,             pset);
  old_dw->get(pdisp,          lb->pDispLabel,          pset);
  old_dw->get(pmass,          lb->pMassLabel,          pset);

  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;
    double X = px[idx].x()-pdisp[idx].x();
    double Y = px[idx].y()-pdisp[idx].y();

    double R = sqrt(X*X + Y*Y);
    double theta = atan2(Y,X);
 
    double alpha = A*sin(c*M_PI*time)*(1. - 32.*pow((R-1),2.) + 256.*pow((R-1),4.));
    double beta = theta + alpha;      

    double br=-R*pow(A,2)*pow(cos(c*M_PI*time),2)*pow(c,2)*pow(M_PI,2)*pow((4*R-3),4)*pow((4*R-5),4)+
      (2/rho0)*2048*pow((4*R-5),2)*pow((R-1),2)*pow((4*R-3),2)*pow(A,2)*pow(sin(c*M_PI*time),2)*R*mu;
    double bt=-A*sin(c*M_PI*time)*(R*pow(c,2)*pow(M_PI,2)*pow((4*R-3),2)*pow((4*R-5),2)+
                                   (1/rho0)*64*mu*(96*pow(R,3)-240*pow(R,2)+188*R-45));
    double bz=0.;

    double bx = cos(beta)*br - sin(beta)*bt;
    double by = sin(beta)*br + cos(beta)*bt;

    pExtForce[idx] = (pmass[idx]*Vector(bx,by,bz));
  }
}

//=====================================================================================
// Axis aligned
//=====================================================================================
void 
MMS::initAxisAligned(const MPMFlags* flags,
                     particleIndex pidx,
                     const Point& p,
                     const Vector& dxcc,
                     const Matrix3& size,
                     ParticleVariable<double>& pvolume,
                     ParticleVariable<double>& pmass,
                     ParticleVariable<Point>& position,
                     ParticleVariable<Vector>& pvelocity,
                     ParticleVariable<Vector>& pdisp,
                     ParticleVariable<Matrix3>& pSize)
{
  /* Vector dx = patch->dCell();             // you need to normalize the variable A by the 
     double normalization = dx.length();    // cell spacing so the Linear interpolation will work
     double A=1e-2 * normalization; */
  double A=0.05;
  double mu = 3846.15;
  double bulk = 8333.33;
  double E = 9.*bulk*mu/(3.*bulk+mu);
  double rho0=1.0;
  double c = sqrt(E/rho0);
  double U = 0.0;
  double V = A*sin(M_PI*p.y())*sin(M_PI*(2./3.));
  double W = A*sin(M_PI*p.z())*sin(M_PI*(4./3.));
  double Fxx=1;
  double Fyy=1+A*M_PI*cos(M_PI*p.y())*sin(2./3.*M_PI);
  double Fzz=1+A*M_PI*cos(M_PI*p.z())*sin(4./3.*M_PI);
  double J = Fxx*Fyy*Fzz;
  double u = A*(c*M_PI)*sin(M_PI*p.x())*cos(0.);
  double v = A*(c*M_PI)*sin(M_PI*p.y())*cos(M_PI*(2./3.));
  double w = A*(c*M_PI)*sin(M_PI*p.z())*cos(M_PI*(4./3.));

  Vector disp(U,V,W);

  position[pidx] = p+disp;
  if (flags->d_axisymmetric) {
    // assume unit radian extent in the circumferential direction
    pvolume[pidx]  = p.x()*(size(0,0)*size(1,1)-size(0,1)*size(1,0))*dxcc.x()*dxcc.y();
  } else {
    // standard voxel volume
    pvolume[pidx]  = J*size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  }

  pSize[pidx] = Matrix3(size(0,0)*Fxx,0.,0.,
                        0.,size(1,1)*Fyy,0.,
                        0.,0.,size(2,2)*Fzz);

  pvelocity[pidx]    = Vector(u,v,w);
  pmass[pidx]        = rho0*pvolume[pidx]/J;
  pdisp[pidx]        = disp;
}

void 
MMS::extForceAxisAligned(const MPMFlags* flags,
                         const MPMLabel* lb,
                         const double& time,
                         ParticleSubset* pset,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw,
                         ParticleVariable<Vector>& pExtForce)
{
  // std::cout << "Entered the AA loop " << std::endl;
  double mu = 3846.15;
  double bulk = 8333.33;
  double E = 9.*bulk*mu/(3.*bulk+mu);
  double lam = (3.*bulk-2.*mu)/3.;
  double rho0=1.0;
  double c = sqrt(E/rho0);
  double A=0.05;

  constParticleVariable<Point>  px;
  constParticleVariable<Vector> pdisp;
  constParticleVariable<double> pmass;
  old_dw->get(px,             lb->pXLabel,             pset);
  old_dw->get(pdisp,          lb->pDispLabel,          pset);
  old_dw->get(pmass,          lb->pMassLabel,          pset);

  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;
    double X = px[idx].x()-pdisp[idx].x();
    double Y = px[idx].y()-pdisp[idx].y();
    double Z = px[idx].z()-pdisp[idx].z();

    double U = A*sin(M_PI*X)*sin(M_PI*c*time);
    double V = A*sin(M_PI*Y)*sin(M_PI*((2./3.)+c*time));
    double W = A*sin(M_PI*Z)*sin(M_PI*((4./3.)+c*time));

    double Fxx=1.+A*M_PI*cos(M_PI*X)*sin(c*M_PI*time);
    double Fyy=1.+A*M_PI*cos(M_PI*Y)*sin(M_PI*((2./3.)+c*time));
    double Fzz=1.+A*M_PI*cos(M_PI*Z)*sin(M_PI*((4./3.)+c*time));

    double K = log(Fxx*Fyy*Fzz);

    double bx=(M_PI*M_PI)*U*(mu/rho0-c*c-(lam*(K-1.)-mu)/(rho0*Fxx*Fxx));
    double by=(M_PI*M_PI)*V*(mu/rho0-c*c-(lam*(K-1.)-mu)/(rho0*Fyy*Fyy));
    double bz=(M_PI*M_PI)*W*(mu/rho0-c*c-(lam*(K-1.)-mu)/(rho0*Fzz*Fzz));

    pExtForce[idx] = pmass[idx]*Vector(bx,by,bz);
  }

}

//=====================================================================================
// Expanding ring
//=====================================================================================
void 
MMS::initExpandingRing(const MPMFlags* flags,
                       particleIndex pidx,
                       const Point& p,
                       const Vector& dxcc,
                       const Matrix3& size,
                       ParticleVariable<double>& pvolume,
                       ParticleVariable<double>& pmass,
                       ParticleVariable<Point>& position,
                       ParticleVariable<Vector>& pvelocity,
                       ParticleVariable<Vector>& pdisp,
                       ParticleVariable<Matrix3>& pSize)
{
  // std::cout << "Entered the ER loop " << std::endl;
  double A = 0.1;
  double E = 1e7;
  double rho = 1000;;
  double C = sqrt(E/rho);
  double ri = 0.4;
  double ro = 0.6;
  double c1 = (-6.*ri)/(ro*(ro - 3.*ri));
  double c2 = (3.*(ro + ri))/(pow(ro,2)*(ro - 3.*ri));
  double c3 = -2./(pow(ro,2)*(ro - 3.*ri));
  double R = sqrt(p.x()*p.x() + p.y()*p.y());

  double u = A*C*M_PI*(c3*pow(R,2.) + c2*R + c1)*p.x();
  double v = A*C*M_PI*(c3*pow(R,2.) + c2*R + c1)*p.y();
  double w = 0.0;
  position[pidx] = p;
  if (flags->d_axisymmetric) {
    // assume unit radian extent in the circumferential direction
    pvolume[pidx]  = p.x()*(size(0,0)*size(1,1)-size(0,1)*size(1,0))*dxcc.x()*dxcc.y();
  } else {
    // standard voxel volume
    pvolume[pidx]  = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  }
  pSize[pidx]        = size;
  pvelocity[pidx]    = Vector(u,v,w);
  pmass[pidx]        = rho*pvolume[pidx];
  pdisp[pidx]        = Vector(0.,0.,0.);
}

void 
MMS::extForceExpandingRing(const MPMFlags* flags,
                           const MPMLabel* lb,
                           const double& time,
                           ParticleSubset* pset,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw,
                           ParticleVariable<Vector>& pExtForce)
{
  double A = 0.1;
  double E = 1e7;
  double mu = E/2.0;
  double rho  = 1000.0;
  double C = sqrt(E/rho);
  double ri = 0.4;
  double ro = 0.6;
  double rb = (ri + ro)/2;
  double c1 = (-6.*ri)/(ro*(ro - 3.*ri));
  double c2 = (3.*(ro + ri))/(pow(ro,2)*(ro - 3.*ri));
  double c3 = -2./(pow(ro,2)*(ro - 3.*ri)); 

  constParticleVariable<Point>  px;
  constParticleVariable<Vector> pdisp;
  constParticleVariable<double> pmass;
  old_dw->get(px,             lb->pXLabel,             pset);
  old_dw->get(pdisp,          lb->pDispLabel,          pset);
  old_dw->get(pmass,          lb->pMassLabel,          pset);

  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;
    double X = px[idx].x()-pdisp[idx].x();
    double Y = px[idx].y()-pdisp[idx].y();
    // std::cout << "time :: " << time << std::endl;
    double t = time;

    double pi_c_t = C*M_PI*t/2.0;
    double pict_rb = pi_c_t/rb;
    double sinct = sin(pict_rb);
    double sinctSq = sinct*sinct;
    double XSq = X*X;
    double YSq = Y*Y;
    double lenSq = XSq + YSq;
    double len = std::sqrt(lenSq);
    double ASq = A*A;
    double c1c2 = c1*c2;
    double c1c3 = c1*c3;
    double c2c3 = c2*c3;
    double c1Sq = c1*c1;
    double c2Sq = c2*c2;
    double c3Sq = c3*c3;
    double Asinct = A*sinct;
    double XAsinct = X*Asinct;
    double YAsinct = Y*Asinct;
    double XYAsinct = X*YAsinct;
    double Xc3Xc2 = X*c3*2.0+X*c2*1.0/len;
    double Yc3Yc2 = Y*c3*2.0+Y*c2*1.0/len;
    double c1c3c2 = c1+c3*lenSq+c2*len;
    double c2pc3 = c2+c3*len*2.0;
    double term1 = Asinct*c1c3c2+XAsinct*Xc3Xc2+1.0;
    double term2 = Asinct*c1c3c2+YAsinct*Yc3Yc2+1.0;
    double term3 = A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0;
    double term4 = c3*Asinct*4.0+ASq*c2Sq*sinctSq*2.0+ASq*c3Sq*sinctSq*lenSq*3.0+ASq*c1c3*sinctSq*4.0;

    // std::cout << "t :: " << t << std::endl;
    // Don't get overwhelmed by looking at the body force. Its insanely large.
    double bx = -((mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(XSq*c3*YAsinct*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-XSq*YAsinct*c2pc3*term3*2.0-YAsinct*c2pc3*lenSq*term3+YAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)+XYAsinct*(X*term4*2.0+ASq*X*c3Sq*sinctSq*lenSq*6.0)*c2pc3*len-XSq*c3*YAsinct*len*term3*2.0+XSq*YAsinct*c2pc3*1.0/len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)-(A*A*A)*XSq*Y*c2c3*pow(sinct,3.0)*c2pc3*lenSq*1.0E1))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))+
                  (mu*(pow(term2,2.0)+ASq*YSq*sinctSq*pow(Xc3Xc2,2.0)-1.0)*(YSq*c3*XAsinct*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-A*X*YSq*sinct*c2pc3*term3*2.0-XAsinct*c2pc3*lenSq*term3+XAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)+XYAsinct*(Y*term4*2.0+ASq*Y*c3Sq*sinctSq*lenSq*6.0)*c2pc3*len-YSq*c3*XAsinct*len*term3*2.0+A*X*YSq*sinct*c2pc3*1.0/len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)-(A*A*A)*X*YSq*c2c3*pow(sinct,3.0)*c2pc3*lenSq*1.0E1))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))-
                  (mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*(Asinct*Yc3Yc2*term2+YAsinct*Xc3Xc2*(Asinct*Xc3Xc2*2.0+XAsinct*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0)))+YAsinct*term1*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0))+XAsinct*(Asinct*Xc3Xc2-A*X*YSq*c2*sinct*1.0/pow(lenSq,3.0/2.0))*Yc3Yc2-A*XSq*Y*c2*sinct*1.0/pow(lenSq,3.0/2.0)*term2))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))-
                  (mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*((Asinct*Yc3Yc2*2.0+YAsinct*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0)))*term2*2.0+ASq*Y*sinctSq*pow(Xc3Xc2,2.0)*2.0-ASq*X*(YSq*Y)*c2*sinctSq*Xc3Xc2*1.0/pow(lenSq,3.0/2.0)*2.0))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))+
                  (mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(Y*1.0/len+A*Y*c2*sinct*4.0+A*Y*c1*sinct*1.0/len+A*Y*c3*sinct*len*6.0+A*(YSq*Y)*c3*sinct*1.0/len*3.0+XSq*c3*YAsinct*1.0/len))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*(Asinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*(Asinct*Yc3Yc2*2.0+YAsinct*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0)))+XAsinct*term2*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0))+YAsinct*(Asinct*Yc3Yc2-A*XSq*Y*c2*sinct*1.0/pow(lenSq,3.0/2.0))*Xc3Xc2-A*X*YSq*c2*sinct*1.0/pow(lenSq,3.0/2.0)*term1)*(len+A*XSq*c2*sinct+A*YSq*c2*sinct*2.0+c1*Asinct*len+A*XSq*c3*sinct*len+A*YSq*c3*sinct*len*3.0))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*(pow(term1,2.0)+ASq*XSq*sinctSq*pow(Yc3Yc2,2.0)-1.0)*(X*1.0/len+A*X*c2*sinct*2.0+A*X*c1*sinct*1.0/len+A*X*c3*sinct*len*2.0+A*(XSq*X)*c3*sinct*1.0/len+YSq*c3*XAsinct*1.0/len*3.0))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*((Asinct*Xc3Xc2*2.0+XAsinct*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0)))*term1*2.0+ASq*X*sinctSq*pow(Yc3Yc2,2.0)*2.0-ASq*(XSq*X)*Y*c2*sinctSq*Yc3Yc2*1.0/pow(lenSq,3.0/2.0)*2.0)*(len+A*XSq*c2*sinct+A*YSq*c2*sinct*2.0+c1*Asinct*len+A*XSq*c3*sinct*len+A*YSq*c3*sinct*len*3.0))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)-
                  mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(len+A*XSq*c2*sinct+A*YSq*c2*sinct*2.0+c1*Asinct*len+A*XSq*c3*sinct*len+A*YSq*c3*sinct*len*3.0)*(Y*1.0/len+A*Y*c2*sinct*6.0+ASq*Y*c1c2*sinctSq*6.0+A*Y*c1*sinct*1.0/len*2.0+A*Y*c3*sinct*len*1.2E1+ASq*Y*c1Sq*sinctSq*1.0/len+ASq*Y*c2Sq*sinctSq*len*6.0+ASq*Y*c3Sq*sinctSq*pow(lenSq,3.0/2.0)*1.5E1+ASq*Y*c1c3*sinctSq*len*1.2E1+ASq*Y*c2c3*sinctSq*lenSq*2.0E1)*1.0/pow(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0,2.0)-
                  mu*(pow(term1,2.0)+ASq*XSq*sinctSq*pow(Yc3Yc2,2.0)-1.0)*(len+A*XSq*c2*sinct+A*YSq*c2*sinct*2.0+c1*Asinct*len+A*XSq*c3*sinct*len+A*YSq*c3*sinct*len*3.0)*(X*1.0/len+A*X*c2*sinct*6.0+ASq*X*c1c2*sinctSq*6.0+A*X*c1*sinct*1.0/len*2.0+A*X*c3*sinct*len*1.2E1+ASq*X*c1Sq*sinctSq*1.0/len+ASq*X*c2Sq*sinctSq*len*6.0+ASq*X*c3Sq*sinctSq*pow(lenSq,3.0/2.0)*1.5E1+ASq*X*c1c3*sinctSq*len*1.2E1+ASq*X*c2c3*sinctSq*lenSq*2.0E1)*1.0/pow(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0,2.0)-
                  mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*1.0/pow(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0),2.0)*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(X*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0)*2.0-X*lenSq*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)*4.0+(X*term4*2.0+ASq*X*c3Sq*sinctSq*lenSq*6.0)*lenSq*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-ASq*X*c2c3*sinctSq*pow(lenSq,2.0)*term3*2.0E1)-
                  mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*1.0/pow(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0),2.0)*(pow(term2,2.0)+ASq*YSq*sinctSq*pow(Xc3Xc2,2.0)-1.0)*(Y*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0)*2.0-Y*lenSq*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)*4.0+(Y*term4*2.0+ASq*Y*c3Sq*sinctSq*lenSq*6.0)*lenSq*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-ASq*Y*c2c3*sinctSq*pow(lenSq,2.0)*term3*2.0E1))/rho-A*(C*C)*(M_PI*M_PI)*X*1.0/(rb*rb)*sinct*c1c3c2*(1.0/4.0);
    double by = -((mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(YSq*c3*XAsinct*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-A*X*YSq*sinct*c2pc3*term3*2.0-XAsinct*c2pc3*lenSq*term3+XAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)+XYAsinct*(Y*term4*2.0+ASq*Y*c3Sq*sinctSq*lenSq*6.0)*c2pc3*len-YSq*c3*XAsinct*len*term3*2.0+A*X*YSq*sinct*c2pc3*1.0/len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)-(A*A*A)*X*YSq*c2c3*pow(sinct,3.0)*c2pc3*lenSq*1.0E1))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))+
                  (mu*(pow(term1,2.0)+ASq*XSq*sinctSq*pow(Yc3Yc2,2.0)-1.0)*(XSq*c3*YAsinct*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-XSq*YAsinct*c2pc3*term3*2.0-YAsinct*c2pc3*lenSq*term3+YAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)+XYAsinct*(X*term4*2.0+ASq*X*c3Sq*sinctSq*lenSq*6.0)*c2pc3*len-XSq*c3*YAsinct*len*term3*2.0+XSq*YAsinct*c2pc3*1.0/len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)-(A*A*A)*XSq*Y*c2c3*pow(sinct,3.0)*c2pc3*lenSq*1.0E1))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))-
                  (mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*(Asinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*(Asinct*Yc3Yc2*2.0+YAsinct*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0)))+XAsinct*term2*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0))+YAsinct*(Asinct*Yc3Yc2-A*XSq*Y*c2*sinct*1.0/pow(lenSq,3.0/2.0))*Xc3Xc2-A*X*YSq*c2*sinct*1.0/pow(lenSq,3.0/2.0)*term1))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))-
                  (mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*((Asinct*Xc3Xc2*2.0+XAsinct*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0)))*term1*2.0+ASq*X*sinctSq*pow(Yc3Yc2,2.0)*2.0-ASq*(XSq*X)*Y*c2*sinctSq*Yc3Yc2*1.0/pow(lenSq,3.0/2.0)*2.0))/(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0))+
                  (mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(X*1.0/len+A*X*c2*sinct*4.0+A*X*c1*sinct*1.0/len+A*X*c3*sinct*len*6.0+A*(XSq*X)*c3*sinct*1.0/len*3.0+YSq*c3*XAsinct*1.0/len))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*(Asinct*Yc3Yc2*term2+YAsinct*Xc3Xc2*(Asinct*Xc3Xc2*2.0+XAsinct*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0)))+YAsinct*term1*(c3*2.0+c2*1.0/len-XSq*c2*1.0/pow(lenSq,3.0/2.0))+XAsinct*(Asinct*Xc3Xc2-A*X*YSq*c2*sinct*1.0/pow(lenSq,3.0/2.0))*Yc3Yc2-A*XSq*Y*c2*sinct*1.0/pow(lenSq,3.0/2.0)*term2)*(len+A*XSq*c2*sinct*2.0+A*YSq*c2*sinct+c1*Asinct*len+A*XSq*c3*sinct*len*3.0+A*YSq*c3*sinct*len))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*(pow(term2,2.0)+ASq*YSq*sinctSq*pow(Xc3Xc2,2.0)-1.0)*(Y*1.0/len+A*Y*c2*sinct*2.0+A*Y*c1*sinct*1.0/len+A*Y*c3*sinct*len*2.0+A*(YSq*Y)*c3*sinct*1.0/len+XSq*c3*YAsinct*1.0/len*3.0))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)+
                  (mu*((Asinct*Yc3Yc2*2.0+YAsinct*(c3*2.0+c2*1.0/len-YSq*c2*1.0/pow(lenSq,3.0/2.0)))*term2*2.0+ASq*Y*sinctSq*pow(Xc3Xc2,2.0)*2.0-ASq*X*(YSq*Y)*c2*sinctSq*Xc3Xc2*1.0/pow(lenSq,3.0/2.0)*2.0)*(len+A*XSq*c2*sinct*2.0+A*YSq*c2*sinct+c1*Asinct*len+A*XSq*c3*sinct*len*3.0+A*YSq*c3*sinct*len))/(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0)-
                  mu*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(len+A*XSq*c2*sinct*2.0+A*YSq*c2*sinct+c1*Asinct*len+A*XSq*c3*sinct*len*3.0+A*YSq*c3*sinct*len)*(X*1.0/len+A*X*c2*sinct*6.0+ASq*X*c1c2*sinctSq*6.0+A*X*c1*sinct*1.0/len*2.0+A*X*c3*sinct*len*1.2E1+ASq*X*c1Sq*sinctSq*1.0/len+ASq*X*c2Sq*sinctSq*len*6.0+ASq*X*c3Sq*sinctSq*pow(lenSq,3.0/2.0)*1.5E1+ASq*X*c1c3*sinctSq*len*1.2E1+ASq*X*c2c3*sinctSq*lenSq*2.0E1)*1.0/pow(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0,2.0)-
                  mu*(pow(term2,2.0)+ASq*YSq*sinctSq*pow(Xc3Xc2,2.0)-1.0)*(len+A*XSq*c2*sinct*2.0+A*YSq*c2*sinct+c1*Asinct*len+A*XSq*c3*sinct*len*3.0+A*YSq*c3*sinct*len)*(Y*1.0/len+A*Y*c2*sinct*6.0+ASq*Y*c1c2*sinctSq*6.0+A*Y*c1*sinct*1.0/len*2.0+A*Y*c3*sinct*len*1.2E1+ASq*Y*c1Sq*sinctSq*1.0/len+ASq*Y*c2Sq*sinctSq*len*6.0+ASq*Y*c3Sq*sinctSq*pow(lenSq,3.0/2.0)*1.5E1+ASq*Y*c1c3*sinctSq*len*1.2E1+ASq*Y*c2c3*sinctSq*lenSq*2.0E1)*1.0/pow(len+A*c2*sinct*lenSq*3.0+c1*Asinct*len*2.0+c3*Asinct*pow(lenSq,3.0/2.0)*4.0+ASq*c1Sq*sinctSq*len+ASq*c2Sq*sinctSq*pow(lenSq,3.0/2.0)*2.0+ASq*c3Sq*sinctSq*pow(lenSq,5.0/2.0)*3.0+ASq*c1c3*sinctSq*pow(lenSq,3.0/2.0)*4.0+ASq*c2c3*sinctSq*pow(lenSq,2.0)*5.0+ASq*c1c2*sinctSq*lenSq*3.0,2.0)-
                  mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*1.0/pow(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0),2.0)*(YAsinct*Xc3Xc2*term1+XAsinct*Yc3Yc2*term2)*(Y*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0)*2.0-Y*lenSq*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)*4.0+(Y*term4*2.0+ASq*Y*c3Sq*sinctSq*lenSq*6.0)*lenSq*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-ASq*Y*c2c3*sinctSq*pow(lenSq,2.0)*term3*2.0E1)-
                  mu*(XYAsinct*c2pc3*lenSq*term3-XYAsinct*c2pc3*len*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0))*1.0/pow(pow(lenSq,2.0)*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)-lenSq*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0),2.0)*(pow(term1,2.0)+ASq*XSq*sinctSq*pow(Yc3Yc2,2.0)-1.0)*(X*pow(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0,2.0)*2.0-X*lenSq*pow(A*c2*sinct*3.0+ASq*c1c2*sinctSq*3.0+ASq*c2c3*sinctSq*lenSq*5.0,2.0)*4.0+(X*term4*2.0+ASq*X*c3Sq*sinctSq*lenSq*6.0)*lenSq*(lenSq*term4+c1*Asinct*2.0+ASq*c1Sq*sinctSq+1.0)*2.0-ASq*X*c2c3*sinctSq*pow(lenSq,2.0)*term3*2.0E1))/rho-A*(C*C)*(M_PI*M_PI)*Y*1.0/(rb*rb)*sinct*c1c3c2*(1.0/4.0);
    double bz = 0.0;
    pExtForce[idx] = pmass[idx]*Vector(bx,by,bz);
  }
}

//=====================================================================================
// Uniaxial strain with harmonic displacement
//=====================================================================================
void 
MMS::initUniaxialStrainHarmonic(const MPMFlags* flags,
                                particleIndex pidx,
                                const Point& p,
                                const Vector& dxcc,
                                const Matrix3& size,
                                ParticleVariable<double>& pvolume,
                                ParticleVariable<double>& pmass,
                                ParticleVariable<Point>& position,
                                ParticleVariable<Vector>& pvelocity,
                                ParticleVariable<Vector>& pdisp,
                                ParticleVariable<Matrix3>& pSize)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  double rho0 = 1700.0;
  double kappa = 6.0e7;
  double mu = 1.0e8;
  double lambda = kappa - (mu*2.0)/3.0;
  double cp = std::sqrt((lambda + 2.0*mu)/rho0);

  // Hardcoded amplitude (0.01 m) and frequency (10000 rad/s)
  double alpha = 0.01;
  double omega = 1000.0;

  // Compute initial velocity and displacement
  double X = p.x();
  double u0 = alpha*std::cos(omega*X/cp);
  double v0 = omega*alpha*std::sin(omega*X/cp);

  // Initialize particle variables
  pvolume[pidx]   = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  pmass[pidx]     = rho0*pvolume[pidx];
  pvelocity[pidx] = Vector(v0, 0.0, 0.0);
  pdisp[pidx]     = Vector(u0, 0.0, 0.0);
  position[pidx]  = p + pdisp[pidx];
  pSize[pidx]     = size;
}

//=====================================================================================
// Body force : Uniaxial strain with non-zero initial stress
//=====================================================================================
void 
MMS::bodyForceUniaxialStrainHarmonic(const MPMLabel* lb,
                                     const double& time,
                                     ParticleSubset* pset,
                                     DataWarehouse* old_dw,
                                     ParticleVariable<Vector>& pBodyForce)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  double rho0 = 1700.0;
  double kappa = 6.0e7;
  double mu = 1.0e8;
  double cp = std::sqrt((kappa + 4.0/3.0*mu)/rho0);

  // Hardcoded amplitude (0.01 m) and frequency (10000 rad/s)
  double alpha = 0.01;
  double omega = 1000.0;

  double omegaSq = omega*omega;
  double omegat = omega*time;
  double omega_by_cp = omega/cp;

  // Get the required particle data
  constParticleVariable<double> pMass;
  constParticleVariable<Point>  pPos;
  constParticleVariable<Vector> pDisp;
  old_dw->get(pMass,          lb->pMassLabel,          pset);
  old_dw->get(pPos,           lb->pXLabel,             pset);
  old_dw->get(pDisp,          lb->pDispLabel,          pset);
  //constParticleVariable<int>    pLoadCurveID;
  //old_dw->get(pLoadCurveID,   lb->pLoadCurveIDLabel,   pset);

  // Loop through particles
  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Skip particles which already have external forces assigned to them
    //int loadCurveID = pLoadCurveID[idx]-1;
    //if (loadCurveID < 0) {

    // Compute body force
    double X = pPos[idx].x() - pDisp[idx].x();
    double omegaX_by_cp = omega_by_cp*X;

    double Urt = alpha*std::cos(omegat - omegaX_by_cp);
    double Uit = -alpha*std::sin(omegat - omegaX_by_cp);

    double bodyForce = (omegaSq*cp*Urt/(cp - omega*Uit) - omegaSq*Urt)*pMass[idx];

    /*
    if (std::abs(pBodyForce[idx].x()) > 0.0) {
      std::cout << "X = " << X << " U = " << pDisp[idx].x() 
                << " m = " << pMass[idx] 
                << " BodyForce = " << bodyForce 
                << " fext_old = " << pBodyForce[idx].x() 
                << " fext_new = " << pBodyForce[idx].x() + bodyForce << std::endl;
    } else {
      if (idx == 100) {
        std::cout << "X = " << X << " U = " << pDisp[idx].x() 
                  << " Urt = " << Urt << " Uit = " << Uit
                  << " m = " << pMass[idx] << " cp = " << cp
                  << " BodyForce = " << bodyForce 
                  << " fext_old = " << pBodyForce[idx].x() 
                  << " fext_new = " << pBodyForce[idx].x() + bodyForce << std::endl;
      }
    }
    */

    pBodyForce[idx] += Vector(bodyForce, 0.0, 0.0);
    //}
  }
}

void 
MMS::initUniaxialStrainHomogeneousLinear(const MPMFlags* flags,
                                         particleIndex pidx,
                                         const Point& p,
                                         const Vector& dxcc,
                                         const Matrix3& size,
                                         ParticleVariable<double>& pvolume,
                                         ParticleVariable<double>& pmass,
                                         ParticleVariable<Point>& position,
                                         ParticleVariable<Vector>& pvelocity,
                                         ParticleVariable<Vector>& pdisp,
                                         ParticleVariable<Matrix3>& pSize)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  double rho0 = 1700.0;
  //double kappa = 6.0e7;
  //double mu = 1.0e8;

  // Amplitude
  double alpha0 = 0.01;

  // Compute initial velocity and displacement
  double X = p.x();
  double u0 = 0.0;
  double v0 = alpha0*X;

  // Initialize particle variables
  pvolume[pidx]   = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  pmass[pidx]     = rho0*pvolume[pidx];
  pvelocity[pidx] = Vector(v0, 0.0, 0.0);
  pdisp[pidx]     = Vector(u0, 0.0, 0.0);
  position[pidx]  = p + pdisp[pidx];
  pSize[pidx]     = size;
}

void 
MMS::bodyForceUniaxialStrainHomogeneousLinear(const MPMLabel* lb,
                                              const double& time,
                                              ParticleSubset* pset,
                                              DataWarehouse* old_dw,
                                              ParticleVariable<Vector>& pBodyForce)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  // double rho0 = 1700.0;
  // double kappa = 6.0e7;
  // double mu = 1.0e8;

  // Hardcoded amplitude (0.01 m)
  // double alpha0 = 0.01;

  // Loop through particles
  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Compute body force
    pBodyForce[idx] += Vector(0.0, 0.0, 0.0);
  }
}

void 
MMS::initUniaxialStrainHomogeneousQuadratic(const MPMFlags* flags,
                                            particleIndex pidx,
                                            const Point& p,
                                            const Vector& dxcc,
                                            const Matrix3& size,
                                            ParticleVariable<double>& pvolume,
                                            ParticleVariable<double>& pmass,
                                            ParticleVariable<Point>& position,
                                            ParticleVariable<Vector>& pvelocity,
                                            ParticleVariable<Vector>& pdisp,
                                            ParticleVariable<Matrix3>& pSize)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  double rho0 = 1700.0;
  //double kappa = 6.0e7;
  //double mu = 1.0e8;

  // Amplitude
  double alpha0 = 0.01;

  // Compute initial velocity and displacement
  double t = 0.0;
  double X = p.x();
  double u0 = alpha0*X*t*t;
  double v0 = 2.0*alpha0*X*t;

  // Initialize particle variables
  pvolume[pidx]   = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  pmass[pidx]     = rho0*pvolume[pidx];
  pvelocity[pidx] = Vector(v0, 0.0, 0.0);
  pdisp[pidx]     = Vector(u0, 0.0, 0.0);
  position[pidx]  = p + pdisp[pidx];
  pSize[pidx]     = size;
}

void 
MMS::bodyForceUniaxialStrainHomogeneousQuadratic(const MPMLabel* lb,
                                                 const double& time,
                                                 ParticleSubset* pset,
                                                 DataWarehouse* old_dw,
                                                 ParticleVariable<Vector>& pBodyForce)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  // double rho0 = 1700.0;
  // double kappa = 6.0e7;
  // double mu = 1.0e8;

  // Hardcoded amplitude (0.01 m)
  double alpha0 = 0.01;

  // Get the required particle data
  constParticleVariable<double> pMass;
  constParticleVariable<Point>  pPos;
  constParticleVariable<Vector> pDisp;
  old_dw->get(pMass,          lb->pMassLabel,          pset);
  old_dw->get(pPos,           lb->pXLabel,             pset);
  old_dw->get(pDisp,          lb->pDispLabel,          pset);
  //constParticleVariable<int>    pLoadCurveID;
  //old_dw->get(pLoadCurveID,   lb->pLoadCurveIDLabel,   pset);

  // Loop through particles
  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Skip particles which already have external forces assigned to them
    //int loadCurveID = pLoadCurveID[idx]-1;
    //if (loadCurveID < 0) {
    double X = pPos[idx].x() - pDisp[idx].x();

    // Compute body force
    double bodyForce = 2.0*alpha0*X*pMass[idx];
    /*
    if (idx == 100) {
      std::cout << "X = " << X << " U = " << pDisp[idx].x() 
                << " m = " << pMass[idx] 
                << " BodyForce = " << bodyForce 
                << " fext_old = " << pBodyForce[idx].x() 
                << " fext_new = " << pBodyForce[idx].x() + bodyForce << std::endl;
    }
    */
    pBodyForce[idx] += Vector(bodyForce, 0.0, 0.0);
    //}
  }
}
