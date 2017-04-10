#ifndef ROOT6_H
#define ROOT6_H

#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>

namespace dem {
bool root6(REAL coef1[], REAL coef2[], Vec& v);
bool root6(REAL coef1[], REAL coef2[], Vec& v, REAL radius, std::size_t partID1,
           std::size_t partID2);
}
#endif
