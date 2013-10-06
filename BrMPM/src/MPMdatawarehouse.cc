#include <MPMdatawarehouse.h>
#include <MPMsaveutil.h>






using namespace Matiti;

 MPMdatawarehouse::MPMdatawarehouse()
             : d_id(0), d_time(new MPMTime())

 ~MPMdatawarehouse::MPMdatawarehouse()



 void
 MPMdatawarehouse::saveData(double dt, MaterialSPArray matlist)
 {
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.saveData(d_out.outputFileCount(), *matlist));
        }
   incrementTime(dt);
   d_id+=1;
 }

        
                 
