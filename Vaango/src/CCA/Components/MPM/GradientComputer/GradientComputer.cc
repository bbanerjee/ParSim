#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>

using namespace Uintah;

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

GradientComputer::GradientComputer(MPMFlags* Mflag)
{
  flag = Mflag;
}

GradientComputer::GradientComputer(const GradientComputer* cm)
{
  flag = cm->flag;
}

GradientComputer::~GradientComputer()
{
}


/*! Calculate gradient of a vector field for 8 noded interpolation */
void 
GradientComputer::computeGrad(Matrix3& grad,
                              vector<IntVector>& ni,
                              vector<Vector>& d_S,
                              const double* oodx, 
                              constNCVariable<Vector>& gVec)
{
  // Compute gradient matrix
  grad.set(0.0);
  for(int k = 0; k < flag->d_8or27; k++) {
    const Vector& vec = gVec[ni[k]];
    for (int j = 0; j<3; j++){
      double fac = d_S[k][j]*oodx[j];
      for (int i = 0; i<3; i++) {
        grad(i,j) += vec[i]*fac;
      }
    }
  }
}

/*! Calculate gradient of vector field for 8 noded interpolation, B matrix
    for Kmat and B matrix for Kgeo */
void 
GradientComputer::computeGradAndBmats(Matrix3& grad,
                                      vector<IntVector>& ni,
                                      vector<Vector>& d_S,
                                      const double* oodx, 
                                      constNCVariable<Vector>& gVec,
                                      const Array3<int>& l2g,
                                      double B[6][24],
                                      double Bnl[3][24],
                                      int* dof)
{
  int l2g_node_num = -1;
      
  computeGrad(grad,ni,d_S,oodx,gVec);
      
  for (int k = 0; k < 8; k++) {
    B[0][3*k] = d_S[k][0]*oodx[0];
    B[3][3*k] = d_S[k][1]*oodx[1];
    B[5][3*k] = d_S[k][2]*oodx[2];
    B[1][3*k] = 0.;
    B[2][3*k] = 0.;
    B[4][3*k] = 0.;
        
    B[1][3*k+1] = d_S[k][1]*oodx[1];
    B[3][3*k+1] = d_S[k][0]*oodx[0];
    B[4][3*k+1] = d_S[k][2]*oodx[2];
    B[0][3*k+1] = 0.;
    B[2][3*k+1] = 0.;
    B[5][3*k+1] = 0.;
        
    B[2][3*k+2] = d_S[k][2]*oodx[2];
    B[4][3*k+2] = d_S[k][1]*oodx[1];
    B[5][3*k+2] = d_S[k][0]*oodx[0];
    B[0][3*k+2] = 0.;
    B[1][3*k+2] = 0.;
    B[3][3*k+2] = 0.;
        
    Bnl[0][3*k] = d_S[k][0]*oodx[0];
    Bnl[1][3*k] = 0.;
    Bnl[2][3*k] = 0.;
    Bnl[0][3*k+1] = 0.;
    Bnl[1][3*k+1] = d_S[k][1]*oodx[1];
    Bnl[2][3*k+1] = 0.;
    Bnl[0][3*k+2] = 0.;
    Bnl[1][3*k+2] = 0.;
    Bnl[2][3*k+2] = d_S[k][2]*oodx[2];
        
    // Need to loop over the neighboring patches l2g to get the right
    // dof number.
    l2g_node_num = l2g[ni[k]];
    dof[3*k]  =l2g_node_num;
    dof[3*k+1]=l2g_node_num+1;
    dof[3*k+2]=l2g_node_num+2;
  }
}

/*! Calculate, for 8 noded interpolation, B matrix
    for Kmat and B matrix for Kgeo */
void 
GradientComputer::computeBmats(vector<IntVector>& ni,
                               vector<Vector>& d_S,
                               const double* oodx, 
                               const Array3<int>& l2g,
                               double B[6][24],
                               double Bnl[3][24],
                               int* dof)
{
  int l2g_node_num = -1;
      
  for (int k = 0; k < 8; k++) {
    B[0][3*k] = d_S[k][0]*oodx[0];
    B[3][3*k] = d_S[k][1]*oodx[1];
    B[5][3*k] = d_S[k][2]*oodx[2];
    B[1][3*k] = 0.;
    B[2][3*k] = 0.;
    B[4][3*k] = 0.;
        
    B[1][3*k+1] = d_S[k][1]*oodx[1];
    B[3][3*k+1] = d_S[k][0]*oodx[0];
    B[4][3*k+1] = d_S[k][2]*oodx[2];
    B[0][3*k+1] = 0.;
    B[2][3*k+1] = 0.;
    B[5][3*k+1] = 0.;
        
    B[2][3*k+2] = d_S[k][2]*oodx[2];
    B[4][3*k+2] = d_S[k][1]*oodx[1];
    B[5][3*k+2] = d_S[k][0]*oodx[0];
    B[0][3*k+2] = 0.;
    B[1][3*k+2] = 0.;
    B[3][3*k+2] = 0.;
        
    Bnl[0][3*k] = d_S[k][0]*oodx[0];
    Bnl[1][3*k] = 0.;
    Bnl[2][3*k] = 0.;
    Bnl[0][3*k+1] = 0.;
    Bnl[1][3*k+1] = d_S[k][1]*oodx[1];
    Bnl[2][3*k+1] = 0.;
    Bnl[0][3*k+2] = 0.;
    Bnl[1][3*k+2] = 0.;
    Bnl[2][3*k+2] = d_S[k][2]*oodx[2];
        
    // Need to loop over the neighboring patches l2g to get the right
    // dof number.
    l2g_node_num = l2g[ni[k]];
    dof[3*k]  =l2g_node_num;
    dof[3*k+1]=l2g_node_num+1;
    dof[3*k+2]=l2g_node_num+2;
  }
}

