#ifndef MATITI_MATERIAL_UP_H
#define MATITI_MATERIAL_UP_H

#include <memory>

namespace Matiti {
  
  class Material;
  typedef std::unique_ptr<Material> MaterialUP;
}

#endif
