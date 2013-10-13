#ifndef __MPMSHAPEFUNCTION_H__
#define __MPMSHAPEFUNCTION_H__

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace MPM
{
  class MPMShapeFunction
  {
  
  public: 

    enum class ShapeType {
       GIMP = 0,
       Quad = 1,
       Linear = 2,
       Cubic = 3
    };
     
    MPMShapeFunction(); 
    ~MPMShapeFunction(); 


    void initialise(const Uintah::ProblemSpecP& ps);

    const ShapeType& ShapeType() const {return d_shape;}
    inline int shapeSize() const {return d_shape_size;}
    inline int numberGhost() const {return d_ghost;}


  private:

    ShapeType d_shape;
    int d_shape_size;
    int d_ghost;

 }; // end class

} // end namespace


#endif
