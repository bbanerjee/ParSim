/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#include <CCA/Components/Schedulers/DependencyException.h>
#include <sstream>

namespace Uintah {

DependencyException::DependencyException(const Task* task,
                                         const VarLabel* label,
                                         int matlIndex,
                                         const Patch* patch,
                                         string has,
                                         string needs,
                                         [[maybe_unused]] const char* file,
                                         [[maybe_unused]] int line)
  : Exception(false)
  , d_task(task)
  , d_label(label)
  , d_mat_index(matlIndex)
  , d_patch(patch)
{
  d_msg = makeMessage(d_task, d_label, d_mat_index, d_patch, has, needs);

#ifdef EXCEPTIONS_CRASH
  std::cout << "A DependencyException exception was thrown.\n";
  std::cout << file << ":" << line << "\n";
  std::cout << d_msg << "\n";
#endif
}

string
DependencyException::makeMessage(const Task* task,
                                 const VarLabel* label,
                                 int matlIndex,
                                 const Patch* patch,
                                 string has,
                                 string needs)
{
  std::ostringstream str;
  str << "Task Dependency Error: (" << has << ") has no corresponding (";
  str << needs << ") for " << label->getName();
  if (patch) {
    str << " on patch " << patch->getID();
    str << ", Level-" << patch->getLevel()->getIndex();
  }

  str << ", for material " << matlIndex;

  if (task != 0) {
    str << " in task " << task->getName();
  }
  str << ".";
  return str.str();
}

DependencyException::DependencyException(const DependencyException& copy)
  : Exception(copy)
  , d_task(copy.d_task)
  , d_label(copy.d_label)
  , d_mat_index(copy.d_mat_index)
  , d_patch(copy.d_patch)
  , d_msg(copy.d_msg)
{
}

const char*
DependencyException::message() const
{
  return d_msg.c_str();
}

const char*
DependencyException::type() const
{
  return "Uintah::Exceptions::DependencyException";
}

} // namespace Uintah
