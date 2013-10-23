#ifndef __MPMSHAPEFUNCTIONFACTORY_H__
#define __MPMSHAPEFUNCTIONFACTORY_H__

#include <ShapeFunctions/MPMShapeFunctionP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace BrMPM
{
  class MPMShapeFunctionFactory
  {
  
  public: 

    MPMShapeFunctionFactory();
    ~MPMShapeFunctionFactory();

    static MPMShapeFunctionP create(const Uintah::ProblemSpecP& ps);

 }; // end class

} // end namespace


#endif
