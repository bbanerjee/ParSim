/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/Core/MPMUtils.h>

namespace Uintah::Util {

auto
face_norm(Patch::FaceType f) -> Vector
{
  switch (f) {
    case Patch::xminus:
      return Vector(-1, 0, 0);
    case Patch::xplus:
      return Vector(1, 0, 0);
    case Patch::yminus:
      return Vector(0, -1, 0);
    case Patch::yplus:
      return Vector(0, 1, 0);
    case Patch::zminus:
      return Vector(0, 0, -1);
    case Patch::zplus:
      return Vector(0, 0, 1);
    default:
      return Vector(0, 0, 0); // oops !
  }
}

auto
getPatchLoHiNodes(const Patch* patch, int n8or27) -> std::pair<IntVector, IntVector>
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

} // end namespace Uintah

