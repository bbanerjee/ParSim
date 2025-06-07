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

#include <StandAlone/Utils/MaterialParticleData.h>
#include <StandAlone/Utils/MaterialParticleVarData.h>

#include <StandAlone/Utils/compare_uda_options.h>
#include <StandAlone/Utils/compare_uda_utils.h>

#include <sstream>

namespace Vaango {
namespace Utils {
namespace CompareUda {

void
MaterialParticleData::createPatchMap()
{
  ASSERT(particleIDs_ != 0); // should check for this before this point
  particleIDs_->createPatchMap();

  std::map<std::string, MaterialParticleVarData>::iterator varIter =
    vars_.begin();
  for (; varIter != vars_.end(); varIter++) {
    (*varIter).second.setParticleIDData(particleIDs_);
  }
}

void
MaterialParticleData::compare(MaterialParticleData& data2,
                              double time1,
                              double time2,
                              double abs_tolerance,
                              double rel_tolerance)
{
  if (vars_.size() == 0) {
    return; // nothing to compare -- all good
  }

  // map particle id's to their patches
  createPatchMap(); // also calls setParticleIDData
  data2.createPatchMap();

  sort();
  data2.sort();

  if (!particleIDs_->compare(*data2.particleIDs_,
                             matl_,
                             time1,
                             time2,
                             abs_tolerance,
                             rel_tolerance)) {
    std::ostringstream warn;
    warn << "    ParticleIDs do not match\n";
    Vaango::Utils::Options::abort_uncomparable(warn);
  }

  std::map<std::string, MaterialParticleVarData>::iterator varIter =
    vars_.begin();
  std::map<std::string, MaterialParticleVarData>::iterator varIter2 =
    data2.vars_.begin();

  for (; (varIter != vars_.end()) && (varIter2 != data2.vars_.end());
       varIter++, varIter2++) {

    // should catch this earlier -- vars/materials do not match
    ASSERT((*varIter).first == (*varIter2).first);

    if ((*varIter).first == "p.particleID") {
      continue; // already compared
    }

    (*varIter).second.compare(
      (*varIter2).second, matl_, time1, time2, abs_tolerance, rel_tolerance);
  }
  // should catch this earlier -- vars/materials do not match
  ASSERT((varIter == vars_.end()) && (varIter2 == data2.vars_.end()));
}

struct Im_Index : public std::pair<Uintah::long64, Uintah::particleIndex>
{
  Im_Index(Uintah::long64 l, Uintah::particleIndex i)
    : std::pair<Uintah::long64, Uintah::particleIndex>(l, i)
  {
  }

  bool
  operator<(Im_Index id2)
  {
    return first < id2.first;
  }
};

void
MaterialParticleData::sort()
{
  // should have made this check earlier -- particleIDs not output
  ASSERT(particleIDs_->getParticleVars().size() != 0);

  std::vector<Im_Index> idIndices;
  Uintah::particleIndex base = 0;

  for (unsigned int i = 0; i < particleIDs_->getParticleVars().size(); i++) {
    Uintah::ParticleVariable<Uintah::long64>* pIDs =
      dynamic_cast<Uintah::ParticleVariable<Uintah::long64>*>(
        particleIDs_->getParticleVars()[i]);

    if (pIDs == 0) {
      std::ostringstream warn;
      warn << "    p.particleID must be a ParticleVariable<long64>\n";
      Vaango::Utils::Options::abort_uncomparable(warn);
    }

    Uintah::long64* pID            = (Uintah::long64*)pIDs->getBasePointer();
    Uintah::ParticleSubset* subset = pIDs->getParticleSubset();

    for (Uintah::ParticleSubset::iterator iter = subset->begin();
         iter != subset->end();
         iter++) {
      idIndices.push_back(Im_Index(*(pID++), base + *iter));
    }
    base = (Uintah::particleIndex)idIndices.size();
  }

  // sort by particle id and find out what happens to the particle indices.
  std::sort(idIndices.begin(), idIndices.end());

  std::vector<Uintah::particleIndex> subsetIndices(idIndices.size());
  for (Uintah::particleIndex i = 0; i < (Uintah::particleIndex)idIndices.size();
       i++) {
    ASSERT(subsetIndices[idIndices[i].second] == 0);
    subsetIndices[idIndices[i].second] = i;
  }

  Uintah::ParticleSubset* subset = scinew Uintah::ParticleSubset(0, matl_, 0);
  subset->expand(subsetIndices.size());
  for (unsigned int i = 0; i < subsetIndices.size(); i++) {
    subset->addParticle(subsetIndices[i]);
  }
  gather(subset);
}

void
MaterialParticleData::gather(Uintah::ParticleSubset* gatherSubset)
{
  std::map<std::string, MaterialParticleVarData>::iterator iter;
  for (iter = vars_.begin(); iter != vars_.end(); iter++) {
    (*iter).second.gather(gatherSubset);
  }
}

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango
