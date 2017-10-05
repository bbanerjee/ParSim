#ifndef ROOT6_OLD_H
#define ROOT6_OLD_H

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

namespace dem {
bool root6_old(REAL coef1[], REAL coef2[], Vec& v, REAL radius,
               std::size_t partID1, std::size_t partID2);
}
#endif
