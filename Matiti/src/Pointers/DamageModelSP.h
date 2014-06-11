#ifndef MATITI_DAMAGE_MODEL_SP_H
#define MATITI_DAMAGE_MODEL_SP_H

#include <memory>

namespace Matiti {
  
  // Forward declaration.  Make sure <DamageModelBase.h> is included before using 
  // DamageModelSP.
  // using stdlib shared_ptr 
  class DamageModelBase;
  typedef std::shared_ptr<DamageModelBase> DamageModelSP;
}

#endif
