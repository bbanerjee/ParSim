#include <cmath>
#include <math.h>

#include <MPMdatawarehouse.h>
#include <MPMsaveutil.h>






using namespace Matiti;

 MPMdatawarehouse::MPMdatawarehouse()
             : d_id(0), d_time(new MPMTime())

 ~MPMdatawarehouse::MPMdatawarehouse()



 void
 MPMdatawarehouse::saveData(double dt, MaterialSPArray& matlist)
 {
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.saveData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id+=1;
 }

 
 void
 MPMdatawarehouse::dumpData(double dt, MaterialSPArray& matlist)
 {
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.dumpData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id+=1;
 }


 bool
 MPMdatawarehouse::checkSave(double dt)
 {
  double dr=d_time.currentTime()/dt;
  double dt0=dt*std::min(dr-std::floor(dr), std::ceil(dr)-dr);
  return dt0<dt/2;
 }

 
 void
 MPMdatawarehouse::init(char lable, int dwi, std::vector val)
 {
  d_id_vec.insert(std::pair<char, std::vector>(lable, val));
 }

 void
 MPMdatawarehouse::append(char lable, int dwi, std::vector val)
 {
  for (auto vec_iter=val.begin(); vec_iter !=val.end(); ++vec_iter) {
            double cur_num = *vec_iter;
            d_id_vec[lable].emplace_back(cur_num);
           }
  }

 
 void
 MPMdatawarehouse::add(char lable, int dwi, std::vector val)
 {
  if (d_id_vec[lable].size() ==0) {
      init(lable, dwi, val);   
   } else {
      append(lable, dwi, val);
   }
 }


 void
 MPMdatawarehouse::zero(char lable, int dwi)
 {
  std::vector zero;
  d_id_vec[lable]=0;  //wrong, it should be modified
 }


 std::vector
 MPMdatawarehouse::get(char lable, int dwi)
 {
  return d_id_vec[lable];
 }


 std::vector < std::vector<double> >
 MPMdatawarehouse::getMult(std::vector<char> lables, int dwi)
 {
   std::vector < std::vector<double> > output;
   for (auto iter = labels.begin(); iter != lables.end(); iter++) {
       char cur_lbl = *iter;
       output.emplace_back(cur_lbl);
   }
   return output;  
 }


 void
 MPMdatawarehouse::addParticles(int dwi, std::vector<MatrixVec>  pointPosition, double pointVolume, double density, int shapeSize)
{
 int const initial = 0.0;
 int const initialOne = 1.0;

 int numberPoints = pointPosition.size();

 std::vector<char> lables = {"pointMomentum", "pointInitialVelocity", "pointInitialPosition", "pointExternalForce",   "pointGradientVelocity", "pointVolumeStress", "pointInternalForce", "pointContactForce", "pointContactMomentum"};

 std::vector<MatrixVec>  pointMomentum, pointInitialVelocity, pointInitialPosition;
 std::vector<MatrixVec>  pointExternalForce, pointInternalForce, pointContactForce;
 std::vector<MatrixVec>  pointContactMomentum, pointMass; 
 
 
        






  
























       
                 
