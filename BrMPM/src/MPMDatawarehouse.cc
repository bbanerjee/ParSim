#include <cmath>
#include <math.h>
#include <stdlib.h>

#include <MPMDatawarehouse.h>
#include <MPMSaveUtil.h>
#include <MPMShapeFunction.h>

using namespace MPM;

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

 
void
MPMDatawarehouse::init(char lable, int dwi, std::vector val)
{
  d_id_vec.insert(std::pair<char, std::vector>(lable, val));
}

void
MPMDatawarehouse::append(char lable, int dwi, std::vector val)
{
  for (auto vec_iter=val.begin(); vec_iter !=val.end(); ++vec_iter) {
            double cur_num = *vec_iter;
            d_id_vec[lable].emplace_back(cur_num);
           }
}

 
void
MPMDatawarehouse::add(char lable, int dwi, std::vector val)
{
  if (d_id_vec[lable].size() ==0) {
      init(lable, dwi, val);   
   } else {
      append(lable, dwi, val);
   }
}


void
MPMDatawarehouse::zero(char lable, int dwi)
{
  std::vector zero;
  d_id_vec[lable]=0;  //wrong, it should be modified
}


std::vector
MPMDatawarehouse::get(char lable, int dwi)
{
  return d_id_vec[lable];
}


std::vector < std::vector<double> >
MPMDatawarehouse::getMult(std::vector<char> lables, int dwi)
{
   std::vector < std::vector<double> > output;
   for (auto iter = labels.begin(); iter != lables.end(); iter++) {
       char cur_lbl = *iter;
       output.emplace_back(cur_lbl);
   }
   return output;  
}


void
MPMDatawarehouse::addParticles(int dwi, MatArrayMatrixVec&  pointsInitialPosition,
                                MatArrayMatrixVec& pointsPosition, 
                                MatArrayMatrixVec& pointsMass, 
                                MatArrayMatrix& pointsGradientVelocity,
                                MatArrayMatrix& pointsStressVelocity, 
                                MatArrayMatrix& pointsDeformationMatrix,
                                MatArrayIntMatrixVecShape& cIndex,
                                MatArrayMatrixVecShape& cWeightFunction,
                                MatArrayMatrixShape& cWeightGradient,
                                std::vector<double>& pointsVolume,
                                double volume, double density)
{
 int Zero = 0;
 double const initialZero = 0.0;
 double const initialOne = 1.0;

 int numberPoints = pointsInitialPosition.size();

 std::vector<char> lables = {"pointMomentum", "pointInitialVelocity", "pointInitialPosition", "pointExternalForce",   "pointGradientVelocity", "pointVolumeStress", "pointInternalForce", "pointContactForce", "pointContactMomentum"};

 initialise(dwi, initialZero, pointsInitialPosition);
 initialise(dwi, initialZero, pointsPosition);
 initialise(dwi, volume*density, pointsMass);
 initialise(dwi, initialZero, pointsGradientVelocity);
 initialise(dwi, initialZero, pointsStressVelocity);
 initialise(dwi, Zero, cIndex);
 initialise(dwi, initialZero, cWeightFunction);
 initialise(dwi, initialZero, cWeightGradient);


 pointsVolume.resize(numberPoints, volume);


 identityMatrix(dwi, initialOne, pointsDeformationMatrix);
} 


void 
MPMDatawarehouse::initialise(int dwi, double initial, MatArrayMarixVec& vec_matrix)
{
  ArrayMatrixVec second_vec_matrix; 
  second_vec_matrix.resize(numberPoints);
  std::vector<MatrixVec>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
      MatrixVec  cur_matrix = *iter;
      cur_matrix(initial);
  }
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}
    
          
 
void 
MPMDatawarehouse::initialise(int dwi, double initial, MatArrayMatrix& vec_matrix)
{
  ArrayMatrix second_vec_matrix;
  second_vec_matrix.resize(numberPoints);
  std::vector<Matrix>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
      Matrix  cur_matrix = *iter;
      cur_matrix(initial);
  }
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}
        

void 
MPMDatawarehouse::initialise(int dwi, int initial, MatArrayIntMatrixVecShape& vec_matrix)
{
  ArrayIntMatrixVecShape second_vec_matrix;
  second_vec_matrix.resize(numberPoints);
  std::vector<IntMatrixVecShape>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
      IntMatrixVecShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}

void 
MPMDatawarehouse::initialise(int dwi, double initial, MatArrayMatrixVecShape& vec_matrix)
{
  ArrayMatrixVecShape second_vec_matrix;
  second_vec_matrix.resize(numberPoints);
  std::vector<MatrixVecShape>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
      MatrixVecShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}

void 
MPMDatawarehouse::initialise(int dwi, double initial, MatArrayMatrixShape& vec_matrix)
{
  ArrayMatrixShape second_vec_matrix;
  second_vec_matrix.resize(numberPoints);
  std::vector<MatrixShape>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
      MatrixShape  cur_matrix = *iter;
      cur_matrix(initial);
  }
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}
 
void 
MPMDatawarehouse::identityMatrix(int dwi, double initial, MatArrayMatrix& vec_matrix)
{
  ArrayMatrix second_vec_matrix;
  second_vec_matrix.resize(numberPoints);
  std::vector<Matrix>::iterator iter;
  for (iter = second_vec_matrix.begin(); iter != second_vec_matrix.end(); iter++) {
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
  vec_matrix = std::make_pair(dwi, second_vec_matrix);
}






  
























       
                 
