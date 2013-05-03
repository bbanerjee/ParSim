#ifndef EMU2DC_MATERIAL_SP_H
#define EMU2DC_MATERIAL_SP_H

#include <memory>

namespace Emu2DC {
  
  // Forward declaration.  Make sure <Node.h> is included before using NodeP.
  // using stdlib shared_ptr 
  class Material;
  typedef std::shared_ptr<Material> MaterialSP;
}

#endif
