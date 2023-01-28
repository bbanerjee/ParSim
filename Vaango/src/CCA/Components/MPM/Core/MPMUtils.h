#ifndef __VAANGO_MPM_UTIL_H__
#define __VAANGO_MPM_UTIL_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Patch.h>

#include <utility>

namespace Uintah {

namespace Util {

std::pair<IntVector, IntVector>
getPatchLoHiNodes(const Patch* patch, int n8or27)
{
  Uintah::IntVector lowIndex(0, 0, 0), highIndex(0, 0, 0);
  if (n8or27 == 8) {
    lowIndex  = patch->getNodeLowIndex();
    highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
  } else if (n8or27 == 27) {
    lowIndex  = patch->getExtraNodeLowIndex();
    highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
  }
  return std::make_pair(lowIndex, highIndex);
}

void
removeDuplicatePatches(Level::selectType& array)
{
  int length = array.size();
  if (length <= 1) {
    return;
  }

  int newLength = 1; // new length of modified array
  int i, j;

  for (i = 1; i < length; i++) {
    for (j = 0; j < newLength; j++) {
      if (array[i] == array[j]) {
        break;
      }
    }
    // if none of the values in array[0..j] == array[i],
    // then copy the current value to a new position in array

    if (j == newLength) {
      array[newLength++] = array[i];
    }
  }
  array.resize(newLength);
}

} // namespace Util
} // end namespace Uintah

#endif // __VAANGO_MPM_UTIL_H__
