#ifndef MATITI_DAMAGEMODEL_UP_H
#define MATITI_DAMAGEMODEL_UP_H

#include <memory>

namespace Matiti {
  
  class DamageModel;
  typedef std::unique_ptr<DamageModel> DamageModelUP;
}

#endif
