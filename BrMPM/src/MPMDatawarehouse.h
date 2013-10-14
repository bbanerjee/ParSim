#ifndef __MATITI_MPMDATAWAREHOUSE__
#define __MATITI_MPMDATAWAREHOUSE__

#include <Geometry/Point3D>
#include <Geometry/Vector3D>
#include <Output.h>
#include <OutputVTK.h>
#include <MPMTime.h>
#include <MPMMatrix.h>
#include <MPMShapeFunction.h>

#include <vector>
#include <iostream>
#include <string>
#include <utility>

namespace MPM {
  
  class MPMDatawarehouse.h {

  public:
    MPMdatawarehouse();
    MPMdatawarehouse(Uintah::ProblemSpecP& ps);
    ~MPMdatawarehouse();

    //void initialise(Uintah::ProblemSpecP& ps);

    void saveData(double dt, MaterialSPArray& matlist);

    void dumpData(double dt, MaterialSPArray& matlist);

    bool checkSave(double dt);

    void init(char label, int dwi, std::vector val);

  private:

    OutputVTK d_out;
    MPMTime d_time;
    MPMSaveUtil d_save;
    MPMShapeFunction d_shapefunction;

    int d_id;
    int const d_dim=3;  //dimension
    int const shapeSize = d_shapefunction.shapeSize();

    typedef std::map<char, std::vector>  vectorIDMap;
    typedef std::vector<std::vector<double>> vectorVector;
    typedef std::vector<std::vector<int>> intVectorVector;

    typedef MPMmatrix<double, 1, d_dim>  MatrixVec;
    typedef std::vector<MatrixVec> ArrayMatrixVec;
    typedef MPMmatrix<double, d_dim, d_dim> Matrix;
    typedef std::vector<Matrix> ArrayMatrix;

    typedef MPMmatrix<int, 1, shapeSize>  IntMatrixVecShape;
    typedef std::vector<IntMatrixVecShape> ArrayIntMatrixVecShape;

    typedef MPMmatrix<double, 1, shapeSize>  MatrixVecShape;
    typedef std::vector<MatrixVecShape> ArrayMatrixVecShape;
    typedef MPMmatrix<double, shapeSize, d_dim> MatrixShape;
    typedef std::vector<MatrixShape> ArrayMatrixShape;

    typedef std::pair<int, ArrayMatrixVec> MatArrayMatrixVec;
    typedef std::pair<int, ArrayMatrix> MatArrayMatrix;
    typedef std::pair<int, ArrayIntMatrixVecShape> MatArrayIntMatrixVecShape;
    typedef std::pair<int, ArrayMatrixVecShape> MatArrayMatrixVecShape;
    typedef std::pair<int, ArrayMatrixShape> MatArrayMatrixShape;

   // typedef std::vector<MPMmatrix<double>> doubleMatrix;
   // typedef std::vector<MPMmatrix<int>>  intMatrix;

   /* std::vector<MatrixVec>  d_pointMomentum, d_pointInitialVelocity, d_pointInitialPosition;
    std::vector<MatrixVec>  d_pointExternalForce, d_pointInternalForce, d_pointContactForce;
    std::vector<MatrixVec>  d_pointContactMomentum, d_pointMass;*/

    vectorIDMap  d_id_vec;

  };
    
} // end namespace

#endif
