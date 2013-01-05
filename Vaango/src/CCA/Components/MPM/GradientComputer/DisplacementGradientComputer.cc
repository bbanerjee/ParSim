#include <CCA/Components/MPM/GradientComputer/DisplacementGradientComputer.h>

using namespace Uintah;

DisplacementGradientComputer::DisplacementGradientComputer(MPMFlags* Mflag) 
  : GradientComputer(Mflag)
{
}

DisplacementGradientComputer::DisplacementGradientComputer(const DisplacementGradientComputer* gc)
  : GradientComputer(gc)
{
}

DisplacementGradientComputer* DisplacementGradientComputer::clone()
{
  return scinew DisplacementGradientComputer(*this);
}

DisplacementGradientComputer::~DisplacementGradientComputer()
{
}

// Actually compute displacement gradient
void 
DisplacementGradientComputer::computeDispGrad(ParticleInterpolator* interp,
                                              const double* oodx,
                                              const Point& px,
                                              const Matrix3& psize,
                                              const Matrix3& pDefGrad_old,
                                              constNCVariable<Vector> gDisp,
                                              Matrix3& dispGrad_new)
{
  // Get the node indices that surround the cell
  vector<IntVector> ni(interp->size());
  vector<Vector> d_S(interp->size());
  interp->findCellAndShapeDerivatives(px, ni, d_S, psize, pDefGrad_old);

  // Compute the gradient
  computeGrad(dispGrad_new, ni, d_S, oodx, gDisp);

  return;
}
