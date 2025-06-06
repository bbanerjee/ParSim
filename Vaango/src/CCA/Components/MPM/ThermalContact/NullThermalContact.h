/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __NullThermalContact__
#define __NullThermalContact__

#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Math/MinMax.h>
#include <cmath>

namespace Uintah {

class ProcessorGroup;
class Patch;
class VarLabel;
class Task;

/**************************************

CLASS
   NullThermalContact

   This version of thermal contact drives the temperatures
   of two materials to the same value at each grid point.

GENERAL INFORMATION

   NullThermalContact.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   NullThermalContact_Model

DESCRIPTION
   Long description...

WARNING

****************************************/

class NullThermalContact : public ThermalContact
{
public:
  // Constructor
  NullThermalContact(ProblemSpecP& ps,
                     const MaterialManagerP& mat_manager,
                     const MPMLabel* labels,
                     const MPMFlags* flags);

  // Destructor
  virtual ~NullThermalContact() = default;

  // Prevent copying/move of this class
  NullThermalContact(const NullThermalContact& con) = delete;
  NullThermalContact(NullThermalContact&& con) = delete;
  NullThermalContact&
  operator=(const NullThermalContact& con) = delete;
  NullThermalContact&
  operator=(NullThermalContact&& con) = delete;

  virtual void
  computeHeatExchange(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw) override;

  virtual void
  initializeThermalContact(const Patch* patch,
                           int vfindex,
                           DataWarehouse* new_dw) override;

  virtual void
  addComputesAndRequires(Task* task,
                         const PatchSet* patches,
                         const MaterialSet* matls) const override;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) override;

};

} // End namespace Uintah

#endif // __NullThermalContact__
