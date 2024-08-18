/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Uintah;

ContactMaterialSpec::ContactMaterialSpec(ProblemSpecP& ps)
{
  if (ps) {
    std::vector<int> materials;
    if (ps->get("materials", materials)) {
      for (const auto material : materials) {
        if (material < 0) {
          throw ProblemSetupException(
            " Invalid material index in contact block", __FILE__, __LINE__);
        }
        this->add(material);
      }
    }
  }
}

void
ContactMaterialSpec::outputProblemSpec(ProblemSpecP& ps)
{
  std::vector<int> matls;
  int i = 0;
  for (const auto material : d_matls) {
    if (material) {
      matls.push_back(i);
    }
    ++i;
  }

  ps->appendElement("materials", matls);
}

void
ContactMaterialSpec::add(unsigned int matlIndex)
{
  // we only add things once at the start, but want
  // quick lookup, so keep logical for each material
  // rather than searching a list every time
  if (d_matls.size() == 0) {
    d_matls.resize(matlIndex + 1);
    for (size_t i = 0; i < matlIndex + 1; i++) {
      d_matls[i] = false;
    }
  }
  if (matlIndex >= d_matls.size()) {
    std::vector<bool> copy(d_matls);
    d_matls.resize(matlIndex + 1);
    for (size_t i = 0; i < copy.size(); i++) {
      d_matls[i] = copy[i];
    }
    for (size_t i = copy.size(); i < matlIndex + 1; i++) {
      d_matls[i] = false;
    }
  }

  d_matls[matlIndex] = true;
}
