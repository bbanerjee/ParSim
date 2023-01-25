/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __COMPOSITE_CONTACT_H__
#define __COMPOSITE_CONTACT_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <list>

namespace Uintah {
using namespace Uintah;

class CompositeContact : public Contact
{
public:
  // Constructor
  CompositeContact(const ProcessorGroup* myworld,
                   const MPMLabel* Mlb,
                   const MPMFlags* MFlag);
  virtual ~CompositeContact();

  void
  outputProblemSpec(ProblemSpecP& ps) override;

  // memory deleted on destruction of composite
  void
  add(std::unique_ptr<Contact> m);

  // how many
  size_t
  size() const
  {
    return d_m.size();
  }

  void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

  void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;

  void
  initFriction(const ProcessorGroup*,
               const PatchSubset*,
               const MaterialSubset* matls,
               DataWarehouse*,
               DataWarehouse* new_dw);

private: // hide
  CompositeContact(const CompositeContact&);
  CompositeContact&
  operator=(const CompositeContact&);

protected: // data
  std::list<std::unique_ptr<Contact>> d_m;
};

} // End namespace Uintah

#endif // __COMPOSITE_CONTACT_H__
