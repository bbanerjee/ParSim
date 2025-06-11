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
#include <CCA/Components/MPM/ReactionDiffusion/SDInterfaceModelFactory.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionInterfaces/CommonIFConcDiff.h>
#include <CCA/Components/MPM/ReactionDiffusion/DiffusionInterfaces/SDInterfaceModel.h>
#include <CCA/Components/MPM/ReactionDiffusion/DiffusionInterfaces/SimpleDiffusionContact.h>
#include <string>

using namespace Uintah;

std::unique_ptr<SDInterfaceModel>
SDInterfaceModelFactory::create(ProblemSpecP& ps,
                                const MaterialManager* ss,
                                const MPMFlags* flags,
                                const MPMLabel* mpm_lb)
{
  ProblemSpecP mpm_ps =
    ps->findBlockWithOutAttribute("MaterialProperties")->findBlock("MPM");
  if (!mpm_ps) {
    throw ProblemSetupException("Cannot MPM material subsection.",
                                __FILE__,
                                __LINE__);
  }

  // Default to a null interface model
  ProblemSpecP child = mpm_ps->findBlock("diffusion_interface");

  // If we don't specify a diffusion interface, assume null.
  std::string diff_interface_type = "null";

  if (child) {
    child->getWithDefault("type", diff_interface_type, "null");
  }

  if (flags->d_integratorType == "implicit") {
    std::string txt = "MPM:  Implicit Scalar Diffusion is not working yet!";
    throw ProblemSetupException(txt, __FILE__, __LINE__);
  }

  if (diff_interface_type == "common") {
    return std::make_unique<CommonIFConcDiff>(child, ss, flags, mpm_lb);
  }

  if (diff_interface_type == "null") {
    return std::make_unique<SDInterfaceModel>(child, ss, flags, mpm_lb);
  }

  if (diff_interface_type == "simple") {
    return std::make_unique<SimpleSDInterface>(child, ss, flags, mpm_lb);
  }

  throw ProblemSetupException("Unknown Scalar Interface Type (" +
                                diff_interface_type + ")",
                              __FILE__,
                              __LINE__);
  return nullptr;
}
