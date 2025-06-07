/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

/*
 *  DependencyException.h
 *
 *  Written by:
 *   Wayne Witzel
 *   Department of Computer Science
 *   University of Utah
 *   May 2002
 *
 */

#ifndef __COMPONENTS_SCHEDULERS_DEPENDENCY_EXCEPTION_H__
#define __COMPONENTS_SCHEDULERS_DEPENDENCY_EXCEPTION_H__

#include <Core/Exceptions/Exception.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <string>

namespace Uintah {

class DependencyException : public Exception
{

public:
  DependencyException(const Task* task,
                      const VarLabel* label,
                      int matlIndex,
                      const Patch* patch,
                      std::string has,
                      std::string needs,
                      const char* file,
                      int line);

  DependencyException(const DependencyException& copy);

  DependencyException& operator=(const DependencyException& copy) = delete;

  virtual ~DependencyException() {}

  static std::string makeMessage(const Task* task,
                                 const VarLabel* label,
                                 int matlIndex,
                                 const Patch* patch,
                                 std::string has,
                                 std::string needs);

  virtual const char* message() const;
  virtual const char* type() const;

private:
  const Task* d_task;
  const VarLabel* d_label;
  int d_mat_index;
  const Patch* d_patch;
  std::string d_msg;
};

} // End namespace Uintah

#endif
