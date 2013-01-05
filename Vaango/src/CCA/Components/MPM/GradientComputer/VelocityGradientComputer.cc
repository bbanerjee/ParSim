#include <CCA/Components/MPM/GradientComputer/VelocityGradientComputer.h>
#include <Core/Exceptions/InvalidValue.h>

using namespace Uintah;

VelocityGradientComputer::VelocityGradientComputer(MPMFlags* Mflag) 
  : GradientComputer(Mflag)
{
}

VelocityGradientComputer::VelocityGradientComputer(const VelocityGradientComputer* gc)
  : GradientComputer(gc)
{
}

VelocityGradientComputer* VelocityGradientComputer::clone()
{
  return scinew VelocityGradientComputer(*this);
}

VelocityGradientComputer::~VelocityGradientComputer()
{
}

// Actually compute velocity gradient
void 
VelocityGradientComputer::computeVelGrad(ParticleInterpolator* interpolator,
                                         const double* oodx,
                                         const short pgFld[],
                                         const Point& px,
                                         const Matrix3& pSize,
                                         const Matrix3& pDefGrad_old,
                                         constNCVariable<Vector>& gVelocity,
                                         constNCVariable<Vector>& GVelocity,
                                         Matrix3& velGrad_new)
{
  if(!flag->d_axisymmetric){

    // Get the node indices that surround the cell
    vector<IntVector> ni(interpolator->size());
    vector<Vector>    d_S(interpolator->size());
    interpolator->findCellAndShapeDerivatives(px, ni, d_S, pSize, pDefGrad_old);

    // Fracture
    if (flag->d_fracture) {
      // Special vel grad for fracture
      computeVelocityGradient(velGrad_new, ni, d_S, oodx, pgFld, gVelocity, GVelocity);
    } else {
      // Standard 3d vel grad computation
      computeGrad(velGrad_new, ni, d_S, oodx, gVelocity);
    }
  } else {  // axi-symmetric kinematics
    // Get the node indices that surround the cell
    vector<IntVector> ni(interpolator->size());
    vector<double>    S(interpolator->size());
    vector<Vector>    d_S(interpolator->size());
    interpolator->findCellAndWeightsAndShapeDerivatives(px, ni, S, d_S,
                                                        pSize,
                                                        pDefGrad_old);
    // x -> r, y -> z, z -> theta
    computeAxiSymVelocityGradient(velGrad_new, ni, d_S, S, oodx, gVelocity, px);
  } // endif (!flag->d_axisymmetric)

  if (isnan(velGrad_new.Norm())) {
    std::cerr << " velGrad = " << velGrad_new << endl;
    throw InvalidValue("**ERROR**: Nan in velocity gradient value", __FILE__, __LINE__);
  }
 
  return;
}

//-------------------------------------------------------------------------
// Protected methods
//-------------------------------------------------------------------------
void 
VelocityGradientComputer::computeAxiSymVelocityGradient(Matrix3& velGrad,
                                             vector<IntVector>& ni,
                                             vector<Vector>& d_S,
                                             vector<double>& S,
                                             const double* oodx,
                                             constNCVariable<Vector>& gVelocity,
                                             const Point& px)
{
  // x -> r, y -> z, z -> theta
  for(int k = 0; k < flag->d_8or27; k++) {
    Vector gvel = gVelocity[ni[k]];
    for (int j = 0; j<2; j++){
      for (int i = 0; i<2; i++) {
        velGrad(i,j)+=gvel[i] * d_S[k][j] * oodx[j];
      }
    }
    velGrad(2,2) += gvel.x()*d_S[k].z();
  }
}

void 
VelocityGradientComputer::computeVelocityGradient(Matrix3& velGrad,
                                        vector<IntVector>& ni,
                                        vector<Vector>& d_S,
                                        const double* oodx, 
                                        const short pgFld[],
                                        constNCVariable<Vector>& gVelocity,
                                        constNCVariable<Vector>& GVelocity)
{
  Vector gvel(0.,0.,0);
  for(int k = 0; k < flag->d_8or27; k++) {
    if(pgFld[k]==1)  gvel = gVelocity[ni[k]];
    if(pgFld[k]==2)  gvel = GVelocity[ni[k]];
    for (int j = 0; j<3; j++){
      double d_SXoodx = d_S[k][j]*oodx[j];
      for (int i = 0; i<3; i++) {
        velGrad(i,j) += gvel[i] * d_SXoodx;
      }
    }
  }
}
    
