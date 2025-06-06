/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Examples_ModelFactory_h
#define Packages_Uintah_CCA_Components_Examples_ModelFactory_h

#include <Core/Grid/MaterialManagerP.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>

/**************************************

CLASS
   ModelFactory

   ModelFactory simulation

GENERAL INFORMATION

   ModelFactory.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   ModelFactory

DESCRIPTION
   Long description...

WARNING

****************************************/

namespace Uintah {

class ProcessorGroup;
class ModelInterface;

class ModelFactory
{

public:
  // this function has a switch for all known components

  static std::vector<std::unique_ptr<ModelInterface>>
  makeModels(const ProcessorGroup* myworld,
             const MaterialManagerP materialManager,
             const ProblemSpecP& orig_or_restart_ps,
             const ProblemSpecP& prob_spec);
};
} // namespace Uintah

#endif
