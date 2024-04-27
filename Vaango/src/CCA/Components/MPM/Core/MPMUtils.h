#ifndef __VAANGO_MPM_UTIL_H__
#define __VAANGO_MPM_UTIL_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Patch.h>

#include <utility>

namespace Uintah {

namespace Util {

extern Vector
face_norm(Patch::FaceType f);

extern std::pair<IntVector, IntVector>
getPatchLoHiNodes(const Patch* patch, int n8or27);

extern void
removeDuplicatePatches(Level::selectType& array);

} // namespace Util
} // end namespace Uintah

#endif // __VAANGO_MPM_UTIL_H__
