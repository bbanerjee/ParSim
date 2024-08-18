/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <CCA/Components/Schedulers/SendState.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Grid/Variables/PSPatchMatlGhostRange.h>
#include <Core/Grid/Variables/ParticleSubset.h>

#include <Core/Parallel/CrowdMonitor.h>
#include <Core/Parallel/Parallel.h>

namespace {

// Tags for each CrowdMonitor
struct send_subsets_tag
{};

using send_subsets_monitor = Uintah::CrowdMonitor<send_subsets_tag>;

}

namespace Uintah {

SendState::~SendState()
{
  {
    send_subsets_monitor sends_write_lock{
      Uintah::CrowdMonitor<send_subsets_tag>::WRITER
    };
    for (auto iter = d_send_subsets.begin(); iter != d_send_subsets.end();
         iter++) {
      if (iter->second->removeReference()) {
        delete iter->second;
      }
    }
  }
}

ParticleSubset*
SendState::find_sendset(int dest,
                        const Patch* patch,
                        int matlIndex,
                        IntVector low,
                        IntVector high,
                        int dwid /* =0 */) const
{
  {
    send_subsets_monitor sends_write_lock{
      Uintah::CrowdMonitor<send_subsets_tag>::READER
    };

    ParticleSubset* ret;
    auto iter = d_send_subsets.find(
      std::make_pair(PSPatchMatlGhostRange(patch, matlIndex, low, high, dwid),
                     dest));

    if (iter == d_send_subsets.end()) {
      ret = nullptr;
    } else {
      ret = iter->second;
    }
    return ret;
  }
}

void
SendState::add_sendset(ParticleSubset* sendset,
                       int dest,
                       const Patch* patch,
                       int matlIndex,
                       IntVector low,
                       IntVector high,
                       int dwid /*=0*/)
{
  {
    auto iter = d_send_subsets.find(
      std::make_pair(PSPatchMatlGhostRange(patch, matlIndex, low, high, dwid),
                     dest));
    if (iter != d_send_subsets.end()) {
      std::cout << "sendSubset Already exists for sendset:" << *sendset
                << " on patch:" << *patch << " matl:" << matlIndex << std::endl;
      SCI_THROW(InternalError("sendSubset already exists", __FILE__, __LINE__));
    }
    d_send_subsets[std::make_pair(
      PSPatchMatlGhostRange(patch, matlIndex, low, high, dwid),
      dest)] = sendset;
    sendset->addReference();
  }
}

void
SendState::reset()
{
  {
    send_subsets_monitor sends_write_lock{
      Uintah::CrowdMonitor<send_subsets_tag>::WRITER
    };

    d_send_subsets.clear();
  }
}

void
SendState::print()
{
  {
    send_subsets_monitor sends_write_lock{
      Uintah::CrowdMonitor<send_subsets_tag>::READER
    };

    for (auto iter = d_send_subsets.begin(); iter != d_send_subsets.end();
         iter++) {
      std::cout << Parallel::getMPIRank() << ' ' << *(iter->second)
                << " src/dest: " << iter->first.second << std::endl;
    }
  }
}

} // namespace Uintah