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
 *  UnknownVariable.h:
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 2000
 *
 */

#ifndef UINTAH_EXCEPTIONS_UNKNOWNVARIABLE_H
#define UINTAH_EXCEPTIONS_UNKNOWNVARIABLE_H

#include <Core/Exceptions/Exception.h>
#include <string>

namespace Uintah {

using Uintah::Exception;

class Level;
class Patch;

class UnknownVariable : public Exception
{
public:
  UnknownVariable(const std::string& varname,
                  int dwid,
                  const Patch* patch,
                  int matlIndex,
                  const std::string& extramsg,
                  const char* file,
                  int line);
  UnknownVariable(const std::string& varname,
                  int dwid,
                  const Level* level,
                  int matlIndex,
                  const std::string& extramsg,
                  const char* file,
                  int line);
  UnknownVariable(const UnknownVariable&);
  virtual ~UnknownVariable();

  virtual const char*
  message() const override;

  virtual const char*
  type() const override;

protected:
private:
  std::string d_msg;
  UnknownVariable&
  operator=(const UnknownVariable&);
};
} // End namespace Uintah

#endif
