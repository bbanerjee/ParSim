#ifndef __MATITI_MPMDATAWAREHOUSE__
#define __MATITI_MPMDATAWAREHOUSE__


#include <Geometry/Point3D>
#include <Geometry/Vector3D>
#include <Output.h>
#include <OutputVTK.h>
#include <MPMTime.h>
#include <MPMmatrix.h>

#include <vector>
#include <iostream>
#include <string>
#include <>

namespace Matiti {
  
  class MPMdatawarehouse.h {

  public:
     MPMdatawarehouse();
     ~MPMdatawarehouse();

  void
 MPMdatawarehouse::saveData(double dt, MaterialSPArray& matlist);

 void
 MPMdatawarehouse::dumpData(double dt, MaterialSPArray& matlist);

 bool
 MPMdatawarehouse::checkSave(double dt);

 void
 MPMdatawarehouse::init(char lable, int dwi, std::vector val);



  private:



    OutputVTK d_out;
    MPMTime d_time;
    MPMsaveutil d_save;


    int d_id;
    int const d_dim=3;  //dimension



    typedef std::map<char, std::vector>  vectorIDMap;
    typedef std::vector<std::vector<double>> vectorVector;
    typedef std::vector<std::vector<int>> intVectorVector;
    typedef MPMmatrix<double, 1, d_dim>  MatrixVec;
    typedef std::vector<MatrixVec> ArrayMatrixVec;
    typedef MPMmatrix<double, d_dim, d_dim> Matrix;
    typedef std::vector<Matrix> ArrayMatrix;


   // typedef std::vector<MPMmatrix<double>> doubleMatrix;
   // typedef std::vector<MPMmatrix<int>>  intMatrix;



   /* std::vector<MatrixVec>  d_pointMomentum, d_pointInitialVelocity, d_pointInitialPosition;
    std::vector<MatrixVec>  d_pointExternalForce, d_pointInternalForce, d_pointContactForce;
    std::vector<MatrixVec>  d_pointContactMomentum, d_pointMass;*/


    vectorIDMap  d_id_vec;
    























#endif
