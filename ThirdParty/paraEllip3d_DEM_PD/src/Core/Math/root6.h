#ifndef ROOT6_H
#define ROOT6_H

#include <Core/Types/realtypes.h>
#include <Core/Math/Vec.h>

namespace dem {
  bool root6(REAL coef1[], REAL coef2[], Vec & v);
}
#endif
