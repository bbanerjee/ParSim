#ifndef __MATITI_MPMDATAWAREHOUSE__
#define __MATITI_MPMDATAWAREHOUSE__

#include <Output.h>
#include <OutputVTK.h>
#include <MPMTime.h>
#include <MPMMatrix.h>
#include <MPMShapeFunction.h>
#include <MPMParticleData.h>
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <map>

#include <boost/variant.hpp>

namespace BrMPM {

class MPMDatawarehouse {

  // Define the variants
  boost::variant<IntegerParticleData, DoubleParticleData, StringParticleData,
                 PointParticleData, VectorParticleData, MatrixParticleData>
             MPMParticleVar;

  public:
    MPMDatawarehouse();
    MPMDatawarehouse(Uintah::ProblemSpecP& ps);
    ~MPMDatawarehouse();

    void saveData(double dt, MaterialSPArray& matlist);

    void dumpData(double dt, MaterialSPArray& matlist);

    bool checkSave(double dt);

    template<typename T>
    void init(char label, int dwi, std::vector<T>& val);

    template<typename T>
    void append(const std::string& label,
                int dwi,
                const std::vector<T>& val);

    template<typename T>
    void add(const std::string& label,
             int dwi,
             const std::vector<T>& val);

    void zero(const std::string& label,
              int dwi);

    template<typename T>
    void get(const std::string& label,
             int dwi,
             std::vector<T>& val);

    template<typename T>
    void getMult(const std::vector<std::string>& labels,
                 int dwi,
                 std::vector<std::vector<T> >& output);

    void addParticles(int dwi,
                      std::vector<Point3D>& pX,
                      std::vector<double>& pVol,
                      std::vector<Vector3D>& pN,
                      std::vector<double>& density,
                      std::vector<int>& shSize);

    void createGrid(int dwi,
                    MPMPatchP& patch);

    void zeroGrid(int dwi);

  private:

    // Datawarehouse data map
    // The string contains the join of the variable name and the 
    // material index (dwi)
    std::map<std::string, MPMParticleVar> d_dw;

  OutputVTK d_out;
  MPMTime d_time;
  MPMSaveUtil d_save;
  MPMShapeFunction d_shapefunction;

  int d_id;
  int const d_dim=3;  //dimension
  int const shapeSize = d_shapefunction.shapeSize();

  typedef std::vector<std::vector<double>> vectorVector;
  typedef std::vector<std::vector<int>> intVectorVector;

  typedef MPMMatrix<double, 1, d_dim>  MatrixVec;
  typedef std::vector<MatrixVec> ArrayMatrixVec;
  typedef MPMMatrix<double, d_dim, d_dim> Matrix;
  typedef std::vector<Matrix> ArrayMatrix;

  typedef MPMMatrix<int, 1, shapeSize>  IntMatrixVecShape;
  typedef std::vector<IntMatrixVecShape> ArrayIntMatrixVecShape;

  typedef MPMMatrix<double, 1, shapeSize>  MatrixVecShape;
  typedef std::vector<MatrixVecShape> ArrayMatrixVecShape;
  typedef MPMMatrix<double, shapeSize, d_dim> MatrixShape;
  typedef std::vector<MatrixShape> ArrayMatrixShape;

    typedef std::pair<int, ArrayMatrixVec> MatArrayMatrixVec;
    typedef std::pair<int, ArrayMatrix> MatArrayMatrix;
    typedef std::pair<int, ArrayIntMatrixVecShape> MatArrayIntMatrixVecShape;
    typedef std::pair<int, ArrayMatrixVecShape> MatArrayMatrixVecShape;
    typedef std::pair<int, ArrayMatrixShape> MatArrayMatrixShape;

  /* std::vector<MatrixVec>  d_pointMomentum, d_pointInitialVelocity, d_pointInitialPosition;
    std::vector<MatrixVec>  d_pointExternalForce, d_pointInternalForce, d_pointContactForce;
    std::vector<MatrixVec>  d_pointContactMomentum, d_pointMass;*/

  vectorIDMap  d_id_vec;

};

}

 // end namespace

#endif
