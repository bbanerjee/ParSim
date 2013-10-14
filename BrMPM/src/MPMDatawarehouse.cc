#include <cmath>
#include <math.h>
#include <stdlib.h>

#include <MPMDatawarehouse.h>
#include <MPMSaveUtil.h>
#include <MPMShapeFunction.h>

using namespace BrMPM;

MPMDatawarehouse::MPMDatawarehouse()
             : d_id(0), d_time(new MPMTime())
{
 /* d_pointMomentum.reserve(1000);
  d_pointInitialVelocity.reserve(1000);
  d_pointInitialPosition.reserve(1000);
  d_pointExternalForce.reserve(1000);
  d_pointInternalForce.reserve(1000);
  d_pointContactForce.reserve(1000);
  d_pointContactMomentum.reserve(1000);
  d_pointMass.reserve(1000); */
}

MPMDatawarehouse::MPMDatawarehouse(Uintah::ProblemSpecP& ps)
{
   d_shapefunction.initialise(ps);
}


MPMDatawarehouse::~MPMDatawarehouse() {}


/* void
 MPMDatawarehouse::initialise(Uintah::ProblemSpecP& ps)
 {
   d_shapefunction.initialise(ps);
 } */


void
MPMDatawarehouse::saveData(double dt, MaterialSPArray& matlist)
{
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.saveData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id += 1;
}

 
void
MPMDatawarehouse::dumpData(double dt, MaterialSPArray& matlist)
{
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.dumpData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id+=1;
}


bool
MPMDatawarehouse::checkSave(double dt)
{
  double dr=d_time.currentTime()/dt;
  double dt0=dt*std::min(dr-std::floor(dr), std::ceil(dr)-dr);
  return dt0<dt/2;
}

 
template<typename T>
inline void MPMDatawarehouse::init(char label, int dwi,
                                   std::vector<T>& val)
{
}

void
MPMDatawarehouse::init(char label, int dwi, std::vector val)
{
  d_id_vec.insert(std::pair<char, std::vector>(label, val));
}

void
MPMDatawarehouse::append(char label, int dwi, std::vector val)
{
  for (auto vec_iter=val.begin(); vec_iter !=val.end(); ++vec_iter) {
            double cur_num = *vec_iter;
            d_id_vec[label].emplace_back(cur_num);
           }
}

 
void
MPMDatawarehouse::add(char label, int dwi, std::vector val)
{
  if (d_id_vec[label].size() ==0) {
      init(label, dwi, val);
   } else {
      append(label, dwi, val);
   }
}


void
MPMDatawarehouse::zero(char label, int dwi)
{
  std::vector zero;
  d_id_vec[label]=0;  //wrong, it should be modified
}


std::vector
MPMDatawarehouse::get(char label, int dwi)
{
  return d_id_vec[label];
}


std::vector < std::vector<double> >
MPMDatawarehouse::getMult(std::vector<char> labels, int dwi)
{
   std::vector < std::vector<double> > output;
   for (auto iter = labels.begin(); iter != labels.end(); iter++) {
       char cur_lbl = *iter;
       output.emplace_back(cur_lbl);
   }
   return output;  
}


void
MPMDatawarehouse::addParticles(int dwi, ArrayMatrixVec&  pointsInitialPosition,
                                ArrayMatrixVec& pointsPosition, 
                                ArrayMatrixVec& pointsMass, 
                                ArrayMatrix& pointsGradientVelocity,
                                ArrayMatrix& pointsStressVelocity, 
                                ArrayMatrix& pointsDeformationMatrix,
                                ArrayIntMatrixVecShape& cIndex,
                                ArrayMatrixVecShape& cWeightFunction,
                                ArrayMatrixShape& cWeightGradient,
                                std::vector<double>& pointsVolume,
                                double volume, double density)
{
 int Zero = 0;
 double const initialZero = 0.0;
 double const initialOne = 1.0;

 int numberPoints = pointsInitialPosition.size();

 std::vector<char> labels = {"pointMomentum", "pointInitialVelocity", "pointInitialPosition", "pointExternalForce",   "pointGradientVelocity", "pointVolumeStress", "pointInternalForce", "pointContactForce", "pointContactMomentum"};

 initialise(initialZero, pointsInitialPosition);
 initialise(initialZero, pointsPosition);
 initialise(volume*density, pointsMass);
 initialise(initialZero, pointsGradientVelocity);
 initialise(initialZero, pointsStressVelocity);
 initialise(Zero, cIndex);
 initialise(initialZero, cWeightFunction);
 initialise(initialZero, cWeightGradient);


 pointsVolume.resize(numberPoints, volume);


 identityMatrix(initialOne, pointsDeformationMatrix);
} 


void 
MPMDatawarehouse::initialise(double initial, ArrayMarixVec& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      MatrixVec  cur_matrix = *iter;
      cur_matrix(initial);
  }
}
    
          
 
void 
MPMDatawarehouse::initialise(double initial, ArrayMatrix& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      Matrix  cur_matrix = *iter;
      cur_matrix(initial);
  }
}
        

void 
MPMDatawarehouse::initialise(int initial,  ArrayIntMatrixVecShape& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      IntMatrixVecShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
}

void 
MPMDatawarehouse::initialise(double initial, ArrayMatrixVecShape& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      MatrixVecShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
}

void 
MPMDatawarehouse::initialise(double initial, ArrayMatrixShape& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      MatrixShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
}
 
void 
MPMDatawarehouse::identityMatrix(double initial, ArrayMatrix& vec_matrix)
{
  vec_matrix.resize(numberPoints);
  for (auto iter = vec_matrix.begin(); iter != vec_matrix.end(); iter++) {
      Matrix  cur_matrix = *iter;
      for (auto mat_iter = cur_matrix.begin(); mat_iter != cur_matrix.end(), mat_iter++) {
          cur_index = *mat_iter;
          div_t divresult;
          divresult = div (cur_index, d_dim);
          int quotion = divresult.quot;
          int remainder = divresult.rem;
          if (quotion == remainder) {
             cur_matrix.set(quotion, remainder, initial);
          }
          else {
             cur_matrix.set(quotion, remainder, 0.0);
          }
       
       } 
                   
  }
}

void MPM::MPMDatawarehouse::zero(const std::string& label, int dwi) {
}

void MPM::MPMDatawarehouse::addParticles(int dwi, std::vector<Point3D>& pX,
    std::vector<double>& pVol, std::vector<Vector3D>& pN,
    std::vector<double>& density, std::vector<int>& shSize) {
}

void MPM::MPMDatawarehouse::createGrid(int dwi, MPMPatchP& patch) {
}

void MPM::MPMDatawarehouse::zeroGrid(int dwi) {
}

//------------------------------------------------------------------------------
// Instantiate templates
template<int>
void init(char label, int dwi, std::vector<int>& val);
template<int>
void append(const std::string& label, int dwi, const std::vector<int>& val);
template<int>
void add(const std::string& label, int dwi, const std::vector<int>& val);
template<int>
void get(const std::string& label, int dwi, std::vector<int>& val);
template<int>
void getMult(const std::vector<std::string>& labels, int dwi,
    std::vector<std::vector<int> >& output);

template<double>
void init(char label, int dwi, std::vector<double>& val);
template<double>
void append(const std::string& label, int dwi, const std::vector<double>& val);
template<double>
void add(const std::string& label, int dwi, const std::vector<double>& val);
template<double>
void get(const std::string& label, int dwi, std::vector<double>& val);
template<double>
void getMult(const std::vector<std::string>& labels, int dwi,
    std::vector<std::vector<double> >& output);

template<std::string>
void init(char label, int dwi, std::vector<std::string>& val);
template<std::string>
void append(const std::string& label, int dwi, const std::vector<std::string>& val);
template<std::string>
void add(const std::string& label, int dwi, const std::vector<std::string>& val);
template<std::string>
void get(const std::string& label, int dwi, std::vector<std::string>& val);
template<std::string>
void getMult(const std::vector<std::string>& labels, int dwi,
    std::vector<std::vector<std::string> >& output);

template<Point3D>
void init(char label, int dwi, std::vector<Point3D>& val);
template<Point3D>
void append(const std::string& label, int dwi, const std::vector<Point3D>& val);
template<Point3D>
void add(const std::string& label, int dwi, const std::vector<Point3D>& val);
template<Point3D>
void get(const std::string& label, int dwi, std::vector<Point3D>& val);
template<Point3D>
void getMult(const std::vector<std::string>& labels, int dwi,
    std::vector<std::vector<Point3D> >& output);

template<Vector3D>
void init(char label, int dwi, std::vector<Vector3D>& val);
template<Vector3D>
void append(const std::string& label, int dwi, const std::vector<Vector3D>& val);
template<Vector3D>
void add(const std::string& label, int dwi, const std::vector<Vector3D>& val);
template<Vector3D>
void get(const std::string& label, int dwi, std::vector<Vector3D>& val);
template<Vector3D>
void getMult(const std::vector<std::string>& labels, int dwi,
    std::vector<std::vector<Vector3D> >& output);

//------------------------------------------------------------------------------






  
























       
                 
