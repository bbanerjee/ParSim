#ifndef __VAANGO_MPM_UTIL_H__
#define __VAANGO_MPM_UTIL_H__

#include <Core/Grid/Patch.h>
#include <Core/Geometry/IntVector.h>

#include <utility>

namespace Uintah {

namespace Util {

std::pair<IntVector, IntVector>
getPatchLoHiNodes(const Patch* patch, int n8or27)
{
   Uintah::IntVector lowIndex(0,0,0), highIndex(0,0,0);
   if (n8or27 == 8) {
     lowIndex = patch->getNodeLowIndex();
     highIndex = patch->getNodeHighIndex()+IntVector(1,1,1);
   } else if (n8or27 == 27) {
     lowIndex = patch->getExtraNodeLowIndex();
     highIndex = patch->getExtraNodeHighIndex()+IntVector(1,1,1);
   }
   return std::make_pair(lowIndex, highIndex);
}

} // end namespace Uintah::Util
} // end namespace Uintah

#endif // __VAANGO_MPM_UTIL_H__
