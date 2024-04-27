/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __VAANGO_STANDALONE_UTILS_COMPARE_UDA_UTILS_H__
#define __VAANGO_STANDALONE_UTILS_COMPARE_UDA_UTILS_H__

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeUtils.h>
#include <Core/Grid/Patch.h>
#include <Core/Parallel/Parallel.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace Vaango {
namespace Utils {
namespace CompareUda {

template<class T>
void
print(std::ostream& out, const T& t);

void
print(std::ostream& out, const Uintah::Vector& t);

void
print(std::ostream& out, const Uintah::Point& t);

void
displayProblemLocation(std::ostream& out,
                       const std::string& var,
                       int matl,
                       const Uintah::Patch* patch,
                       double time);

void
displayProblemLocation(std::ostream& out,
                       const std::string& var,
                       int matl,
                       const Uintah::Patch* patch,
                       const Uintah::Patch* patch2,
                       double time);

bool
compare(double a, double b, double abs_tolerance, double rel_tolerance);

bool
compare(float a, float b, double abs_tolerance, double rel_tolerance);

bool
compare(Uintah::long64 a,
        Uintah::long64 b,
        double /* abs_tolerance */,
        double /* rel_tolerance */);

bool
compare(int a, int b, double /* abs_tolerance */, double /* rel_tolerance */);

bool
compare(Uintah::Vector a,
        Uintah::Vector b,
        double abs_tolerance,
        double rel_tolerance);

bool
compare(Uintah::IntVector a,
        Uintah::IntVector b,
        double abs_tolerance,
        double rel_tolerance);

bool
compare(Uintah::Point a,
        Uintah::Point b,
        double abs_tolerance,
        double rel_tolerance);

bool
compare(Uintah::Stencil7& a,
        Uintah::Stencil7& b,
        double abs_tolerance,
        double rel_tolerance);

bool
compare(const Uintah::Matrix3& a,
        const Uintah::Matrix3& b,
        double abs_tolerance,
        double rel_tolerance);

template<class T>
void
compareParticles(Uintah::DataArchive* da1,
                 Uintah::DataArchive* da2,
                 const std::string& var_name,
                 int matl,
                 const Uintah::Patch* patch1,
                 const Uintah::Patch* patch2,
                 double time,
                 int timestep,
                 double abs_tolerance,
                 double rel_tolerance);

template<class T>
void
comparePerPatch(Uintah::DataArchive* da1,
                Uintah::DataArchive* da2,
                const std::string& var_name,
                int matl,
                const Uintah::Patch* patch1,
                const Uintah::Patch* patch2,
                double time,
                int timestep,
                double abs_tolerance,
                double rel_tolerance);

// map nodes to their owning patch in a level.
void
buildPatchMap(Uintah::LevelP level,
              const std::string& filebase,
              Uintah::Array3<const Uintah::Patch*>& patchMap,
              double time,
              Uintah::Patch::VariableBasis basis);

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_COMPARE_UDA_UTILS_H__
