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

/*
 *  DimensionMismatch.h: Exception to indicate a failed bounds check
 *
 *  Written by:
 *   McKay Davis
 *   Department of Computer Science
 *   University of Utah
 *   Feb 2005
 *
 */

#ifndef Core_Exceptions_GuiException_h
#define Core_Exceptions_GuiException_h

#include <Core/Exceptions/Exception.h>
#include <string>
using std::string;

namespace Uintah {
class GuiException : public Exception
{
public:
  GuiException(string msg)
    : Exception()
    , msg_(msg)
  {
    stacktrace_ = 0;
  }
  virtual ~GuiException(){};
  virtual const char*
  message() const override
  {
    return msg_.c_str();
  }
  virtual const char*
  type() const override
  {
    return "GuiException";
  }

protected:
private:
  string msg_;
  GuiException&
  operator=(const GuiException);
};
} // End namespace Uintah

#endif // Core_Exceptions_GuiException_h
