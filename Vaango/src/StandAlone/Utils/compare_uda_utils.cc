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

#include <StandAlone/Utils/compare_uda_options.h>
#include <StandAlone/Utils/compare_uda_utils.h>

#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Util/ProgressiveWarning.h>

#include <filesystem>
#include <iostream>

namespace Vaango {
namespace Utils {
namespace CompareUda {

// see Vector specialization below
template<class T>
void
print(std::ostream& out, const T& t)
{
  out << t;
}

template void
print(std::ostream& out, const int& t);

template void
print(std::ostream& out, const float& t);

template void
print(std::ostream& out, const double& t);

template void
print(std::ostream& out, const Uintah::Matrix3& t);

template void
print(std::ostream& out, const Uintah::Stencil7& t);

// must override Vector's output in order to use the ostream's precision
void
print(std::ostream& out, const Uintah::Vector& t)
{
  out << "[" << t.x() << ", " << t.y() << ", " << t.z() << "]";
}

// must override Vector's output in order to use the ostream's precision
void
print(std::ostream& out, const Uintah::Point& t)
{
  out << "[" << t.x() << ", " << t.y() << ", " << t.z() << "]";
}

void
displayProblemLocation(std::ostream& out,
                       const std::string& var,
                       int matl,
                       const Uintah::Patch* patch,
                       double time)
{
  out << "Time: " << time << "\n"
      << "Variable: " << var << "\n"
      << "Material: " << matl << "\n";
  if (patch != 0) {
    out << "Patch: " << patch->getID() << "\n";
  }
}

void
displayProblemLocation(std::ostream& out,
                       const std::string& var,
                       int matl,
                       const Uintah::Patch* patch,
                       const Uintah::Patch* patch2,
                       double time)
{
  out << "Time: " << time << " "
      << "Level: " << patch->getLevel()->getIndex() << " "
      << "Patch1: " << patch->getID() << " "
      << "Patch2: " << patch2->getID() << " "
      << "Material: " << matl << " "
      << "Variable: " << var << std::endl;
}

bool
compare(double a, double b, double abs_tolerance, double rel_tolerance)
{
  // Return false only if BOTH absolute and relative comparisons fail.
  if (std::isnan(a) || std::isnan(b)) {
    return false;
  }

  double max_abs = fabs(a);
  if (std::abs(b) > max_abs) {
    max_abs = fabs(b);
  }
  if (std::abs(a - b) > abs_tolerance) {
    if (max_abs > 0 && (std::abs(a - b) / max_abs) > rel_tolerance) {
      return false;
    } else {
      return true;
    }
  } else {
    return true;
  }
}

bool
compare(float a, float b, double abs_tolerance, double rel_tolerance)
{
  if (std::isnan(a) || std::isnan(b)) {
    return false;
  }

  return compare((double)a, (double)b, abs_tolerance, rel_tolerance);
}

bool
compare(Uintah::long64 a,
        Uintah::long64 b,
        double /* abs_tolerance */,
        double /* rel_tolerance */)
{
  if (std::isnan(a) || std::isnan(b)) {
    return false;
  }

  return (a == b); // longs should use an exact comparison
}

bool
compare(int a, int b, double /* abs_tolerance */, double /* rel_tolerance */)
{
  if (std::isnan(a) || std::isnan(b)) {
    return false;
  }

  return (a == b); // int should use an exact comparison
}

bool
compare(Uintah::Vector a,
        Uintah::Vector b,
        double abs_tolerance,
        double rel_tolerance)
{
  if (std::isnan(a.length()) || std::isnan(b.length())) {
    return false;
  }

  return compare(a.x(), b.x(), abs_tolerance, rel_tolerance) &&
         compare(a.y(), b.y(), abs_tolerance, rel_tolerance) &&
         compare(a.z(), b.z(), abs_tolerance, rel_tolerance);
}

bool
compare(Uintah::IntVector a,
        Uintah::IntVector b,
        double abs_tolerance,
        double rel_tolerance)
{
  //  if(std::isnan(a.length()) || std::isnan(b.length())){
  if (std::isnan(a.x() * a.x() + a.y() * a.y() + a.z() * a.z()) ||
      std::isnan(b.x() * b.x() + b.y() * b.y() + b.z() * b.z())) {
    return false;
  }

  return compare(a.x(), b.x(), abs_tolerance, rel_tolerance) &&
         compare(a.y(), b.y(), abs_tolerance, rel_tolerance) &&
         compare(a.z(), b.z(), abs_tolerance, rel_tolerance);
}

bool
compare(Uintah::Point a,
        Uintah::Point b,
        double abs_tolerance,
        double rel_tolerance)
{
  return compare(a.asVector(), b.asVector(), abs_tolerance, rel_tolerance);
}

bool
compare(Uintah::Stencil7& a,
        Uintah::Stencil7& b,
        double abs_tolerance,
        double rel_tolerance)
{
  return compare(a.p, b.p, abs_tolerance, rel_tolerance) &&
         compare(a.n, b.n, abs_tolerance, rel_tolerance) &&
         compare(a.s, b.s, abs_tolerance, rel_tolerance) &&
         compare(a.e, b.e, abs_tolerance, rel_tolerance) &&
         compare(a.w, b.w, abs_tolerance, rel_tolerance) &&
         compare(a.t, b.t, abs_tolerance, rel_tolerance) &&
         compare(a.b, b.b, abs_tolerance, rel_tolerance);
}

bool
compare(const Uintah::Matrix3& a,
        const Uintah::Matrix3& b,
        double abs_tolerance,
        double rel_tolerance)
{
  if (std::isnan(a.Norm()) || std::isnan(b.Norm())) {
    return false;
  }

  // for (int i = 0; i < 3; i++)
  //   for (int j = 0; j < 3; j++)
  //     if (!compare(a(i,j), b(i, j), abs_tolerance, rel_tolerance))
  // Comparing element by element is overly sensitive to code changes
  // The following is a hopefully more informative metric of agreement
  if (!compare(a.Norm(), b.Norm(), abs_tolerance, rel_tolerance)) {
    return false;
  } else {
    return true;
  }
}

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
                 double rel_tolerance)
{
  Uintah::ParticleVariable<T> var1;
  Uintah::ParticleVariable<T> var2;
  da1->query(var1, var_name, matl, patch1, timestep);
  da2->query(var2, var_name, matl, patch2, timestep);

  Uintah::ParticleSubset* pset1 = var1.getParticleSubset();
  Uintah::ParticleSubset* pset2 = var2.getParticleSubset();

  if (pset1->numParticles() != pset2->numParticles()) {
    std::ostringstream warn;
    warn << "Inconsistent number of particles.\n";
    Vaango::Utils::CompareUda::displayProblemLocation(
      warn, var_name, matl, patch1, time);
    warn << "    " << Vaango::Utils::Options::filebase_1() << " has "
         << pset1->numParticles() << " particles.\n";
    warn << "    " << Vaango::Utils::Options::filebase_2() << " has "
         << pset2->numParticles() << " particles.\n";
    Vaango::Utils::Options::abort_uncomparable(warn);
  }

  Uintah::ParticleSubset::iterator iter1 = pset1->begin();
  Uintah::ParticleSubset::iterator iter2 = pset2->begin();

  bool premature_failure{ false };
  for (; iter1 != pset1->end() && iter2 != pset2->end(); iter1++, iter2++) {
    if (!Vaango::Utils::CompareUda::compare(
          var1[*iter1], var2[*iter2], abs_tolerance, rel_tolerance)) {
      std::cerr << "\nValues differ too much.\n";

      Vaango::Utils::CompareUda::displayProblemLocation(
        std::cerr, var_name, matl, patch1, time);

      std::cerr << Vaango::Utils::Options::filebase_1() << ":\n";
      print(std::cerr, var1[*iter1]);

      std::cerr << std::endl << Vaango::Utils::Options::filebase_2() << ":\n";
      print(std::cerr, var2[*iter2]);
      std::cerr << std::endl;

      Vaango::Utils::Options::tolerance_failure();
      if (Vaango::Utils::Options::concise()) {
        premature_failure = true;
        break;
      }
    }
  }

  // this should be true if both sets are the same size
  if (!premature_failure) {
    ASSERT(iter1 == pset1->end() && iter2 == pset2->end());
  }
}

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
                double rel_tolerance)
{
  Uintah::PerPatch<T> var1;
  Uintah::PerPatch<T> var2;
  da1->query(var1, var_name, matl, patch1, timestep);
  da2->query(var2, var_name, matl, patch2, timestep);

  if (!Vaango::Utils::CompareUda::compare(
        var1, var2, abs_tolerance, rel_tolerance)) {

    std::cerr << "\nValues differ too much.\n";

    Vaango::Utils::CompareUda::displayProblemLocation(
      std::cerr, var_name, matl, patch1, time);

    std::cerr << Vaango::Utils::Options::filebase_1() << ":\n";
    print(std::cerr, var1);

    std::cerr << std::endl << Vaango::Utils::Options::filebase_2() << ":\n";
    print(std::cerr, var2);

    std::cerr << std::endl;
    Vaango::Utils::Options::tolerance_failure();
  }
}

// map nodes to their owning patch in a level.
// Nodes are used because I am assuming that whoever owns the node at
// that index also owns the cell, or whatever face at that same index.
// The same doesn't work if you used cells because nodes can go beyond
// cells (when there is no neighbor on the greater side).
void
buildPatchMap(Uintah::LevelP level,
              const std::string& filebase,
              Uintah::Array3<const Uintah::Patch*>& patchMap,
              double time,
              Uintah::Patch::VariableBasis basis)
{
  const Uintah::PatchSet* allPatches = level->allPatches();
  const Uintah::PatchSubset* patches = allPatches->getUnion();
  if (patches->size() == 0) {
    return;
  }

  Uintah::IntVector bl   = Uintah::IntVector(0, 0, 0);
  Uintah::IntVector low  = patches->get(0)->getExtraLowIndex(basis, bl);
  Uintah::IntVector high = patches->get(0)->getExtraHighIndex(basis, bl);

  for (int i = 1; i < patches->size(); i++) {
    low  = Min(low, patches->get(i)->getExtraLowIndex(basis, bl));
    high = Max(high, patches->get(i)->getExtraHighIndex(basis, bl));
  }

  patchMap.resize(low, high);
  patchMap.initialize(0);

  Uintah::Level::const_patch_iterator iter;
  for (iter = level->patchesBegin(); iter != level->patchesEnd(); iter++) {
    const Uintah::Patch* patch = *iter;

    ASSERT(Uintah::Min(patch->getExtraLowIndex(basis, bl), low) == low);
    ASSERT(Uintah::Max(patch->getExtraHighIndex(basis, bl), high) == high);

    patchMap.rewindow(patch->getExtraLowIndex(basis, bl),
                      patch->getExtraHighIndex(basis, bl));

    Uintah::serial_for(patchMap.range(), [&](int i, int j, int k) {
      if (patchMap(i, j, k) != nullptr) {
        static Uintah::ProgressiveWarning pw(
          "Two patches on the same grid overlap", 10);
        if (pw.invoke()) {
          std::cerr << "Patches " << patch->getID() << " and "
                    << (*iter)->getID() << " overlap on the same file at time "
                    << time << " in " << filebase << " at index "
                    << IntVector(i, j, k) << std::endl;
        }
        // abort_uncomparable();

        // in some cases, we can have overlapping patches, where an extra
        // cell/node overlaps an interior cell/node of another patch.  We prefer
        // the interior one.  if there are two overlapping interior ones (nodes
        // or face centers only), they should have the same value.  However,
        // since this patchMap is also used for cell centered variables give
        // priority to the patch that has this index within its interior cell
        // centered variables
        Uintah::IntVector in_low  = patch->getLowIndex(basis);
        Uintah::IntVector in_high = patch->getHighIndex(basis);

        if (i >= in_low.x() && j >= in_low.y() && k >= in_low.z() &&
            i < in_high.x() && j < in_high.y() && k < in_high.z()) {
          patchMap(i, j, k) = patch;
        }
      } else {
        patchMap(i, j, k) = patch;
      }
    });
  }
  patchMap.rewindow(low, high);
}

template void
compareParticles<int>(Uintah::DataArchive*,
                      Uintah::DataArchive*,
                      const std::string&,
                      int,
                      Uintah::Patch const*,
                      Uintah::Patch const*,
                      double,
                      int,
                      double,
                      double);

template void
compareParticles<float>(Uintah::DataArchive*,
                        Uintah::DataArchive*,
                        const std::string&,
                        int,
                        Uintah::Patch const*,
                        Uintah::Patch const*,
                        double,
                        int,
                        double,
                        double);

template void
compareParticles<double>(Uintah::DataArchive*,
                         Uintah::DataArchive*,
                         const std::string&,
                         int,
                         Uintah::Patch const*,
                         Uintah::Patch const*,
                         double,
                         int,
                         double,
                         double);

template void
compareParticles<Uintah::Point>(Uintah::DataArchive*,
                                Uintah::DataArchive*,
                                const std::string&,
                                int,
                                Uintah::Patch const*,
                                Uintah::Patch const*,
                                double,
                                int,
                                double,
                                double);

template void
compareParticles<Uintah::Vector>(Uintah::DataArchive*,
                                 Uintah::DataArchive*,
                                 const std::string&,
                                 int,
                                 Uintah::Patch const*,
                                 Uintah::Patch const*,
                                 double,
                                 int,
                                 double,
                                 double);

template void
compareParticles<Uintah::IntVector>(Uintah::DataArchive*,
                                    Uintah::DataArchive*,
                                    const std::string&,
                                    int,
                                    Uintah::Patch const*,
                                    Uintah::Patch const*,
                                    double,
                                    int,
                                    double,
                                    double);

template void
compareParticles<Uintah::Matrix3>(Uintah::DataArchive*,
                                  Uintah::DataArchive*,
                                  const std::string&,
                                  int,
                                  Uintah::Patch const*,
                                  Uintah::Patch const*,
                                  double,
                                  int,
                                  double,
                                  double);

template void
comparePerPatch<double>(Uintah::DataArchive*,
                        Uintah::DataArchive*,
                        const std::string&,
                        int,
                        Uintah::Patch const*,
                        Uintah::Patch const*,
                        double,
                        int,
                        double,
                        double);

template void
comparePerPatch<float>(Uintah::DataArchive*,
                       Uintah::DataArchive*,
                       const std::string&,
                       int,
                       Uintah::Patch const*,
                       Uintah::Patch const*,
                       double,
                       int,
                       double,
                       double);

template void
comparePerPatch<int>(Uintah::DataArchive*,
                     Uintah::DataArchive*,
                     const std::string&,
                     int,
                     Uintah::Patch const*,
                     Uintah::Patch const*,
                     double,
                     int,
                     double,
                     double);
template void
comparePerPatch<Uintah::Point>(Uintah::DataArchive*,
                               Uintah::DataArchive*,
                               const std::string&,
                               int,
                               Uintah::Patch const*,
                               Uintah::Patch const*,
                               double,
                               int,
                               double,
                               double);
template void
comparePerPatch<Uintah::Vector>(Uintah::DataArchive*,
                                Uintah::DataArchive*,
                                const std::string&,
                                int,
                                Uintah::Patch const*,
                                Uintah::Patch const*,
                                double,
                                int,
                                double,
                                double);
template void
comparePerPatch<Uintah::IntVector>(Uintah::DataArchive*,
                                   Uintah::DataArchive*,
                                   const std::string&,
                                   int,
                                   Uintah::Patch const*,
                                   Uintah::Patch const*,
                                   double,
                                   int,
                                   double,
                                   double);

template void
comparePerPatch<Uintah::Matrix3>(Uintah::DataArchive*,
                                 Uintah::DataArchive*,
                                 const std::string&,
                                 int,
                                 Uintah::Patch const*,
                                 Uintah::Patch const*,
                                 double,
                                 int,
                                 double,
                                 double);

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango
