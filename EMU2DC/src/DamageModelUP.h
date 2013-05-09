#ifndef EMU2DC_DAMAGEMODEL_UP_H
#define EMU2DC_DAMAGEMODEL_UP_H

#include <memory>

namespace Emu2DC {
  
  class DamageModel;
  typedef std::unique_ptr<DamageModel> DamageModelUP;
}

#endif
