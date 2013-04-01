#ifndef MATITI_MATERIAL_H
#define MATITI_MATERIAL_H

#include <Vaango/Core/ProblemSpec/ProblemSpecP.h>
#include <Core/SimulationState.h>

namespace Matiti {

  class Material {

    public:

      Material(Vaango::ProblemSpecP& ps);
      Material();
      
      virtual ~Material();
      
      virtual void registerBondState(SimulationState* ss);

    private:

      Material(const Material &mat);
      Material& operator=(const Material &mat);
   };
} // End namespace Matiti

#endif // MATITI_MATERIAL_H
