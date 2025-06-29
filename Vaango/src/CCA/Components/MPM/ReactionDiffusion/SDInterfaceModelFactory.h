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

#ifndef _SDINTERFACEMODELFACTORY_H_
#define _SDINTERFACEMODELFACTORY_H_

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionInterfaces/SDInterfaceModel.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <memory>
#include <string>

namespace Uintah {

class SDInterface;
class MPMLabel;
class MPMFlags;
class MaterialManager;

class SDInterfaceModelFactory
{
public:
  // this function has a switch for all known mat_types

  static std::unique_ptr<SDInterfaceModel>
  create(ProblemSpecP& ps,
         const MaterialManager* ss,
         const MPMFlags* flags,
         const MPMLabel* mpm_lb);
};
} // End namespace Uintah

#endif /* _CONSTITUTIVEMODELFACTORY_H_ */
