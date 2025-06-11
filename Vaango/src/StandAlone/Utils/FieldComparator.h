/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef __VAANGO_STANDALONE_UTILS_FIELD_COMPARATOR_H__
#define __VAANGO_STANDALONE_UTILS_FIELD_COMPARATOR_H__

#include <Core/Containers/ConsecutiveRangeSet.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/Array3.h>

#include <string>

namespace Vaango {
namespace Utils {
namespace CompareUda {

class FieldComparator
{

public:
  virtual ~FieldComparator(){};
  virtual void
  compareFields(Uintah::DataArchive* da1,
                Uintah::DataArchive* da2,
                const std::string& var,
                Uintah::ConsecutiveRangeSet matls,
                const Uintah::Patch* patch,
                const Uintah::Array3<const Uintah::Patch*>& patch2Map,
                double time,
                int timestep,
                double abs_tolerance,
                double rel_tolerance) = 0;

  static FieldComparator*
  makeFieldComparator(const Uintah::TypeDescription* td,
                      const Uintah::TypeDescription* subtype,
                      const Uintah::Patch* patch,
                      const bool includeExtraCells);
};

//______________________________________________________________________
//
template<class Field, class Iterator>
class SpecificFieldComparator : public FieldComparator
{
public:
  SpecificFieldComparator(Iterator begin)
    : m_begin(begin)
  {
  }
  virtual ~SpecificFieldComparator() {}

  virtual void
  compareFields(Uintah::DataArchive* da1,
                Uintah::DataArchive* da2,
                const std::string& var,
                Uintah::ConsecutiveRangeSet matls,
                const Uintah::Patch* patch,
                const Uintah::Array3<const Uintah::Patch*>& patch2Map,
                double time,
                int timestep,
                double abs_tolerance,
                double rel_tolerance);

private:
  Iterator m_begin;
};

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_FIELD_COMPARATOR_H__
