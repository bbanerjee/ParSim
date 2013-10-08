#ifndef __MATITI_MPMDATAWAREHOUSE__
#define __MATITI_MPMDATAWAREHOUSE__


#include <Geometry/Point3D>
#include <Geometry/Vector3D>
#include <Output.h>
#include <OutputVTK.h>
#include <MPMTime.h>

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



    typedef std::map<int, std::vector>  vectorIDMap;
    vectorIDMap  d_id_vec;
    























#endif
