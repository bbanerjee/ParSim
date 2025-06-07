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

#include <StandAlone/Utils/MaterialParticleVarData.h>

#include <StandAlone/Utils/compare_uda_options.h>
#include <StandAlone/Utils/compare_uda_utils.h>

#include <sstream>

namespace Vaango {
namespace Utils {
namespace CompareUda {

MaterialParticleVarData::~MaterialParticleVarData()
{
  std::vector<Uintah::ParticleVariableBase*>::iterator iter =
    m_particleVars.begin();
  for (; iter != m_particleVars.end(); iter++) {
    delete *iter;
  }

  if (m_name == "p.particleID") {
    delete m_patchMap;
  }
}

void
MaterialParticleVarData::createPatchMap()
{
  ASSERT(m_name == "p.particleID");
  if (m_patchMap) {
    delete m_patchMap;
  }

  m_patchMap = scinew std::map<Uintah::long64, const Uintah::Patch*>();

  for (unsigned int patch = 0; patch < m_particleVars.size(); patch++) {
    Uintah::particleIndex count =
      m_particleVars[patch]->getParticleSubset()->numParticles();

    Uintah::ParticleVariable<Uintah::long64>* particleID =
      dynamic_cast<Uintah::ParticleVariable<Uintah::long64>*>(
        m_particleVars[patch]);

    if (particleID == 0) {
      std::ostringstream warn;
      warn << "    p.particleID must be a ParticleVariable<long64>\n";
      Vaango::Utils::Options::abort_uncomparable(warn);
    }

    for (int i = 0; i < count; i++) {
      (*m_patchMap)[(*particleID)[i]] = m_patches[patch];
    }
  }
}

void
MaterialParticleVarData::add(Uintah::ParticleVariableBase* pvb,
                             const Uintah::Patch* patch)
{
  m_particleVars.push_back(pvb);
  subsets_.push_back(pvb->getParticleSubset());
  m_patches.push_back(patch);
}

void
MaterialParticleVarData::gather(Uintah::ParticleSubset* gatherSubset)
{
  ASSERT(m_particleVars.size() > 0);
  Uintah::ParticleVariableBase* pvb = m_particleVars[0]->clone();
  pvb->gather(gatherSubset, subsets_, m_particleVars, 0);
  m_particleVars.clear();

  subsets_.clear();
  m_patches.clear();
  add(pvb, 0 /* all patches */);
}

bool
MaterialParticleVarData::compare(MaterialParticleVarData& data2,
                                 int matl,
                                 double time1,
                                 double time2,
                                 double abs_tolerance,
                                 double rel_tolerance)
{
  std::cerr << "\tVariable: " << m_name << ", comparing via particle ids"
            << std::endl;
  ASSERT(m_particleVars.size() == 1 && subsets_.size() == 1 &&
         data2.m_particleVars.size() == 1 && data2.subsets_.size() == 1);

  Uintah::ParticleVariableBase* pvb1 = m_particleVars[0];
  Uintah::ParticleVariableBase* pvb2 = data2.m_particleVars[0];

  // type checks should have been made earlier
  ASSERT(pvb1->virtualGetTypeDescription() ==
         pvb2->virtualGetTypeDescription());

  switch (pvb1->virtualGetTypeDescription()->getSubType()->getType()) {
    case Uintah::TypeDescription::Type::double_type:
      return compare(data2,
                     dynamic_cast<Uintah::ParticleVariable<double>*>(pvb1),
                     dynamic_cast<Uintah::ParticleVariable<double>*>(pvb2),
                     matl,
                     time1,
                     time2,
                     abs_tolerance,
                     rel_tolerance);
    case Uintah::TypeDescription::Type::float_type:
      return compare(data2,
                     dynamic_cast<Uintah::ParticleVariable<float>*>(pvb1),
                     dynamic_cast<Uintah::ParticleVariable<float>*>(pvb2),
                     matl,
                     time1,
                     time2,
                     abs_tolerance,
                     rel_tolerance);
    case Uintah::TypeDescription::Type::long64_type:
      return compare(data2,
                     dynamic_cast<Uintah::ParticleVariable<Uintah::long64>*>(pvb1),
                     dynamic_cast<Uintah::ParticleVariable<Uintah::long64>*>(pvb2),
                     matl,
                     time1,
                     time2,
                     abs_tolerance,
                     rel_tolerance);
    case Uintah::TypeDescription::Type::int_type:
      return compare(data2,
                     dynamic_cast<Uintah::ParticleVariable<int>*>(pvb1),
                     dynamic_cast<Uintah::ParticleVariable<int>*>(pvb2),
                     matl,
                     time1,
                     time2,
                     abs_tolerance,
                     rel_tolerance);
    case Uintah::TypeDescription::Type::Point:
      return compare(
        data2,
        dynamic_cast<Uintah::ParticleVariable<Uintah::Point>*>(pvb1),
        dynamic_cast<Uintah::ParticleVariable<Uintah::Point>*>(pvb2),
        matl,
        time1,
        time2,
        abs_tolerance,
        rel_tolerance);
    case Uintah::TypeDescription::Type::Vector:
      return compare(
        data2,
        dynamic_cast<Uintah::ParticleVariable<Uintah::Vector>*>(pvb1),
        dynamic_cast<Uintah::ParticleVariable<Uintah::Vector>*>(pvb2),
        matl,
        time1,
        time2,
        abs_tolerance,
        rel_tolerance);
    case Uintah::TypeDescription::Type::IntVector:
      return compare(
        data2,
        dynamic_cast<Uintah::ParticleVariable<Uintah::IntVector>*>(pvb1),
        dynamic_cast<Uintah::ParticleVariable<Uintah::IntVector>*>(pvb2),
        matl,
        time1,
        time2,
        abs_tolerance,
        rel_tolerance);
    case Uintah::TypeDescription::Type::Matrix3:
      return compare(
        data2,
        dynamic_cast<Uintah::ParticleVariable<Uintah::Matrix3>*>(pvb1),
        dynamic_cast<Uintah::ParticleVariable<Uintah::Matrix3>*>(pvb2),
        matl,
        time1,
        time2,
        abs_tolerance,
        rel_tolerance);
    default:
      std::cerr << "MaterialParticleVarData::gather: ParticleVariable of "
                   "unsupported type: "
                << pvb1->virtualGetTypeDescription()->getName() << '\n';
      Uintah::Parallel::exitAll(-1);
  }
  return 0;
}

template<class T>
bool
MaterialParticleVarData::compare(MaterialParticleVarData& data2,
                                 Uintah::ParticleVariable<T>* value1,
                                 Uintah::ParticleVariable<T>* value2,
                                 int matl,
                                 double time1,
                                 double /*time2*/,
                                 double abs_tolerance,
                                 double rel_tolerance)
{
  bool passes                   = true;
  Uintah::ParticleSubset* pset1 = value1->getParticleSubset();
  Uintah::ParticleSubset* pset2 = value2->getParticleSubset();

  if (pset1->numParticles() != pset2->numParticles()) {
    std::ostringstream warn;
    warn << "Inconsistent number of particles.\n";

    Vaango::Utils::CompareUda::displayProblemLocation(
      warn, m_name, matl, 0, time1);

    warn << "    " << Vaango::Utils::Options::filebase_1() << " has "
         << pset1->numParticles() << " particles.\n";
    warn << "    " << Vaango::Utils::Options::filebase_2() << " has "
         << pset2->numParticles() << " particles.\n";
    Vaango::Utils::Options::abort_uncomparable(warn);
  }

  // Assumes that the particleVariables are in corresponding order --
  // not necessarily by their particle set order.  This is what the
  // sort/gather achieves.
  for (unsigned int i = 0; i < pset1->numParticles(); i++) {
    if (!(Vaango::Utils::CompareUda::compare(
          (*value1)[i], (*value2)[i], abs_tolerance, rel_tolerance))) {
      if (m_name != "p.particleID") {
        ASSERT(getParticleID(i) == data2.getParticleID(i));
      }

      std::cerr << std::setprecision(18) << std::endl;
      std::cerr << "DIFFERENCE on particle id= " << getParticleID(i)
                << std::endl;

      IntVector origin((int)(getParticleID(i) >> 16) & 0xffff,
                       (int)(getParticleID(i) >> 32) & 0xffff,
                       (int)(getParticleID(i) >> 48) & 0xffff);

      std::cerr << "(Originating from " << origin << ")\n";

      const Uintah::Patch* patch1 = getPatch(i);
      const Uintah::Patch* patch2 = data2.getPatch(i);

      Vaango::Utils::CompareUda::displayProblemLocation(
        std::cerr, m_name, matl, patch1, patch2, time1);

      std::cerr << Vaango::Utils::Options::filebase_1() << ":\n"
                << (*value1)[i] << std::endl;
      std::cerr << Vaango::Utils::Options::filebase_2() << ":\n"
                << (*value2)[i] << std::endl;

      Vaango::Utils::Options::tolerance_failure();
      if (Vaango::Utils::Options::concise()) {
        break;
      }

      passes = false;
    }
  }

  return passes;
}

Uintah::long64
MaterialParticleVarData::getParticleID(Uintah::particleIndex index)
{
  ASSERT(m_particleIDData != 0);
  ASSERT(m_particleIDData->m_particleVars.size() == 1);

  Uintah::ParticleVariable<Uintah::long64>* particleIDs =
    dynamic_cast<Uintah::ParticleVariable<Uintah::long64>*>(
      m_particleIDData->m_particleVars[0]);

  ASSERT(particleIDs != 0);

  return (*particleIDs)[index];
}

const Uintah::Patch*
MaterialParticleVarData::getPatch(Uintah::particleIndex index)
{
  ASSERT(m_patchMap != 0);
  return (*m_patchMap)[getParticleID(index)];
}

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango
