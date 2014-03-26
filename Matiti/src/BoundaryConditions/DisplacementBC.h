#ifndef __MATITI_DISPLACEMENT_BC_H__
#define __MATITI_DISPLACEMENT_BC_H__

namespace Matiti {
  
  class DisplacementBC {
  
  public:

    DisplacementBC();
    ~DisplacementBC();

    void initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes);
    void applyBC(const SCIRun::Vector& disp, NodePArray& surfaceNodes); 

  private:

   // prevent copying
   DisplacementBC(const DisplacementBC& bc);
   DisplacementBC& operator=(const DisplacementBC& bc);

  }; // end class

} // end namespace
#endif

