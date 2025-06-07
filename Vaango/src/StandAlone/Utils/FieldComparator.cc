/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <StandAlone/Utils/FieldComparator.h>

#include <StandAlone/Utils/compare_uda_options.h>
#include <StandAlone/Utils/compare_uda_utils.h>

#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Parallel/Parallel.h>

#include <sstream>

namespace Vaango {
namespace Utils {
namespace CompareUda {

FieldComparator*
FieldComparator::makeFieldComparator(const Uintah::TypeDescription* td,
                                     const Uintah::TypeDescription* subtype,
                                     const Uintah::Patch* patch,
                                     const bool includeExtraCells)
{
  switch (td->getType()) {

    case Uintah::TypeDescription::Type::ParticleVariable:
    case Uintah::TypeDescription::Type::PerPatch:
      // Particles and PerPatch handled differently (and previously)
      break;
    case Uintah::TypeDescription::Type::NCVariable: {
      Uintah::NodeIterator iter = patch->getNodeIterator();

      if (includeExtraCells) {
        iter = patch->getExtraNodeIterator();
      }

      switch (subtype->getType()) {
        case Uintah::TypeDescription::Type::double_type:
          return scinew SpecificFieldComparator<Uintah::NCVariable<double>,
                                                Uintah::NodeIterator>(iter);
        case Uintah::TypeDescription::Type::float_type:
          return scinew SpecificFieldComparator<Uintah::NCVariable<float>,
                                                Uintah::NodeIterator>(iter);
        case Uintah::TypeDescription::Type::int_type:
          return scinew SpecificFieldComparator<Uintah::NCVariable<int>,
                                                Uintah::NodeIterator>(iter);
        case Uintah::TypeDescription::Type::Point:
          return scinew
            SpecificFieldComparator<Uintah::NCVariable<Uintah::Point>,
                                    Uintah::NodeIterator>(iter);
        case Uintah::TypeDescription::Type::Vector:
          return scinew
            SpecificFieldComparator<Uintah::NCVariable<Uintah::Vector>,
                                    Uintah::NodeIterator>(iter);
        case Uintah::TypeDescription::Type::Matrix3:
          return scinew
            SpecificFieldComparator<Uintah::NCVariable<Uintah::Matrix3>,
                                    Uintah::NodeIterator>(iter);
        default:
          std::cerr << "FieldComparator::makeFieldComparator: NC Variable of "
                       "unsupported type: "
                    << subtype->getName() << '\n';
          Uintah::Parallel::exitAll(-1);
      }
    }

    case Uintah::TypeDescription::Type::CCVariable: {
      Uintah::CellIterator iter = patch->getCellIterator();

      if (includeExtraCells) {
        iter = patch->getExtraCellIterator();
      }

      switch (subtype->getType()) {
        case Uintah::TypeDescription::Type::double_type:
          return scinew SpecificFieldComparator<Uintah::CCVariable<double>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::float_type:
          return scinew SpecificFieldComparator<Uintah::CCVariable<float>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::int_type:
          return scinew SpecificFieldComparator<Uintah::CCVariable<int>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Point:
          return scinew
            SpecificFieldComparator<Uintah::CCVariable<Uintah::Point>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Vector:
          return scinew
            SpecificFieldComparator<Uintah::CCVariable<Uintah::Vector>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Matrix3:
          return scinew
            SpecificFieldComparator<Uintah::CCVariable<Uintah::Matrix3>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Stencil7:
          return scinew
            SpecificFieldComparator<Uintah::CCVariable<Uintah::Stencil7>,
                                    Uintah::CellIterator>(iter);
        default:
          std::cerr << "FieldComparator::makeFieldComparator: CC Variable of "
                       "unsupported type: "
                    << subtype->getName() << '\n';
          Uintah::Parallel::exitAll(-1);
      }
    }

    case Uintah::TypeDescription::Type::SFCXVariable: {
      Uintah::CellIterator iter = patch->getSFCXIterator();

      if (includeExtraCells) {
        iter = patch->getExtraSFCXIterator();
      }

      switch (subtype->getType()) {
        case Uintah::TypeDescription::Type::double_type:
          return scinew SpecificFieldComparator<Uintah::SFCXVariable<double>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::float_type:
          return scinew SpecificFieldComparator<Uintah::SFCXVariable<float>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::int_type:
          return scinew SpecificFieldComparator<Uintah::SFCXVariable<int>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Point:
          return scinew
            SpecificFieldComparator<Uintah::SFCXVariable<Uintah::Point>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Vector:
          return scinew
            SpecificFieldComparator<Uintah::SFCXVariable<Uintah::Vector>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Matrix3:
          return scinew
            SpecificFieldComparator<Uintah::SFCXVariable<Uintah::Matrix3>,
                                    Uintah::CellIterator>(iter);
        default:
          std::cerr << "FieldComparator::makeFieldComparator: SFCX Variable of "
                       "unsupported type: "
                    << subtype->getName() << '\n';
          Uintah::Parallel::exitAll(-1);
      }
    }

    case Uintah::TypeDescription::Type::SFCYVariable: {
      Uintah::CellIterator iter = patch->getSFCYIterator();

      if (includeExtraCells) {
        iter = patch->getExtraSFCYIterator();
      }

      switch (subtype->getType()) {
        case Uintah::TypeDescription::Type::double_type:
          return scinew SpecificFieldComparator<Uintah::SFCYVariable<double>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::float_type:
          return scinew SpecificFieldComparator<Uintah::SFCYVariable<float>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::int_type:
          return scinew SpecificFieldComparator<Uintah::SFCYVariable<int>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Point:
          return scinew
            SpecificFieldComparator<Uintah::SFCYVariable<Uintah::Point>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Vector:
          return scinew
            SpecificFieldComparator<Uintah::SFCYVariable<Uintah::Vector>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Matrix3:
          return scinew
            SpecificFieldComparator<Uintah::SFCYVariable<Uintah::Matrix3>,
                                    Uintah::CellIterator>(iter);
        default:
          std::cerr << "FieldComparator::makeFieldComparator: SFCY Variable of "
                       "unsupported type: "
                    << subtype->getName() << '\n';
          Uintah::Parallel::exitAll(-1);
      }
    }

    case Uintah::TypeDescription::Type::SFCZVariable: {
      Uintah::CellIterator iter = patch->getSFCZIterator();

      if (includeExtraCells) {
        iter = patch->getExtraSFCZIterator();
      }

      switch (subtype->getType()) {
        case Uintah::TypeDescription::Type::double_type:
          return scinew SpecificFieldComparator<Uintah::SFCZVariable<double>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::float_type:
          return scinew SpecificFieldComparator<Uintah::SFCZVariable<float>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::int_type:
          return scinew SpecificFieldComparator<Uintah::SFCZVariable<int>,
                                                Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Point:
          return scinew
            SpecificFieldComparator<Uintah::SFCZVariable<Uintah::Point>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Vector:
          return scinew
            SpecificFieldComparator<Uintah::SFCZVariable<Uintah::Vector>,
                                    Uintah::CellIterator>(iter);
        case Uintah::TypeDescription::Type::Matrix3:
          return scinew
            SpecificFieldComparator<Uintah::SFCZVariable<Uintah::Matrix3>,
                                    Uintah::CellIterator>(iter);
        default:
          std::cerr << "FieldComparator::makeFieldComparator: SFCZ Variable of "
                       "unsupported type: "
                    << subtype->getName() << '\n';
          Uintah::Parallel::exitAll(-1);
      }
    }

    default:
      std::cerr
        << "FieldComparator::makeFieldComparator: Variable of unsupported "
           "type: "
        << td->getName() << '\n';
      Uintah::Parallel::exitAll(-1);
  }
  return nullptr;
}

template<class Field, class Iterator>
void
SpecificFieldComparator<Field, Iterator>::compareFields(
  Uintah::DataArchive* da1,
  Uintah::DataArchive* da2,
  const std::string& var_name,
  Uintah::ConsecutiveRangeSet matls,
  const Uintah::Patch* patch,
  const Uintah::Array3<const Uintah::Patch*>& patch2Map,
  double time1,
  int timestep,
  double abs_tolerance,
  double rel_tolerance)
{
  Field* field2;
  bool firstMatl = true;

  //__________________________________
  //  Matl loop
  for (Uintah::ConsecutiveRangeSet::iterator matlIter = matls.begin();
       matlIter != matls.end();
       matlIter++) {

    int matl = *matlIter;

    Field field;
    bool found = da1->query(field, var_name, matl, patch, timestep);

    if (!found) {
      std::cout << "Skipping comparison of " << var_name
                << " as it was not found in DataArchive1.\n";
      continue;
    }

    std::map<const Uintah::Patch*, Field*> patch2FieldMap;
    typename std::map<const Uintah::Patch*, Field*>::iterator findIter;

    for (Iterator iter = m_begin; !iter.done(); iter++) {
      const Uintah::Patch* patch2 = patch2Map[*iter];

      findIter = patch2FieldMap.find(patch2);

      if (findIter == patch2FieldMap.end()) {

        if (firstMatl) { // check only needs to be made the first round
          Uintah::ConsecutiveRangeSet matls2 =
            da2->queryMaterials(var_name, patch2, timestep);
          ASSERT(matls == matls2); // check should have been made previously
        }

        field2                 = scinew Field();
        patch2FieldMap[patch2] = field2;
        found = da2->query(*field2, var_name, matl, patch2, timestep);

        if (!found) {
          std::cout << "Skipping comparison of " << var_name
                    << " as it was not found in DataArchive2.\n";
          continue;
        }
      } else {
        field2 = (*findIter).second;
      }

      if (!Vaango::Utils::CompareUda::compare(
            field[*iter], (*field2)[*iter], abs_tolerance, rel_tolerance)) {

        std::cerr << "DIFFERENCE " << *iter << "  ";
        Vaango::Utils::CompareUda::displayProblemLocation(
          std::cerr, var_name, matl, patch, patch2, time1);

        std::cerr << Vaango::Utils::Options::filebase_1() << " (1)\t\t"
                  << Vaango::Utils::Options::filebase_2() << " (2)"
                  << std::endl;
        Vaango::Utils::CompareUda::print(std::cerr, field[*iter]);
        std::cerr << "\t\t";
        Vaango::Utils::CompareUda::print(std::cerr, (*field2)[*iter]);
        std::cerr << std::endl;

        Vaango::Utils::Options::tolerance_failure();
        if (Vaango::Utils::Options::concise()) {
          break; // Exit for() loop as we are only displaying first error per
                 // variable.
        }
      }
    }

    typename std::map<const Uintah::Patch*, Field*>::iterator iter =
      patch2FieldMap.begin();
    for (; iter != patch2FieldMap.end(); iter++) {
      delete (*iter).second;
    }
    firstMatl = false;
  }
}

template void
SpecificFieldComparator<Uintah::SFCZVariable<int>, Uintah::CellIterator>::
  compareFields(Uintah::DataArchive*,
                Uintah::DataArchive*,
                const std::string&,
                Uintah::ConsecutiveRangeSet,
                Uintah::Patch const*,
                Uintah::Array3<Uintah::Patch const*> const&,
                double,
                int,
                double,
                double);

template void
SpecificFieldComparator<Uintah::NCVariable<Uintah::Stencil7>,
                        Uintah::CellIterator>::
  compareFields(Uintah::DataArchive*,
                Uintah::DataArchive*,
                const std::string&,
                Uintah::ConsecutiveRangeSet,
                Uintah::Patch const*,
                Uintah::Array3<Uintah::Patch const*> const&,
                double,
                int,
                double,
                double);

template void
SpecificFieldComparator<Uintah::CCVariable<Uintah::Stencil7>,
                        Uintah::CellIterator>::
  compareFields(Uintah::DataArchive*,
                Uintah::DataArchive*,
                const std::string&,
                Uintah::ConsecutiveRangeSet,
                Uintah::Patch const*,
                Uintah::Array3<Uintah::Patch const*> const&,
                double,
                int,
                double,
                double);

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango
