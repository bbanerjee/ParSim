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
 *  Exception.h: Base class for all SCI Exceptions
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   July 1999
 *
 */

#ifndef Core_Exceptions_Exception_h
#define Core_Exceptions_Exception_h

#if USE_SCI_THROW
#define SCI_THROW(exc) do {Uintah::Exception::sci_throw(exc);throw exc;} while(Uintah::Exception::alwaysFalse())
#else
#define SCI_THROW(exc) throw exc
#endif

#include <string>

#include <sci_defs/error_defs.h>


namespace Uintah {

  ///////////////////////////////////////////////////////////////////////////////////
  // If "wait for debugger" is turned on, then anytime an exception is thrown,
  // the code will 'pause', print out a hostname/PID, and wait until a debugger
  // is attached.  (Defaults to off.)  

  void TURN_ON_WAIT_FOR_DEBUGGER();
  void TURN_OFF_WAIT_FOR_DEBUGGER();
  void WAIT_FOR_DEBUGGER(bool useFlag=false); // Note, if turned off, this call does nothing!

  //
  ///////////////////////////////////////////////////////////////////////////////////

  class  Exception {
  public:
    Exception(bool ignoreWait=false);
    virtual ~Exception();
    virtual const char* message() const=0;
    virtual const char* type() const=0;
    const char* stackTrace() const {
      return stacktrace_;
    }

    static void sci_throw(const Exception& exc);
    static bool alwaysFalse();
  protected:
    const char* stacktrace_;
  private:
    Exception& operator=(const Exception&);
  };

   std::string getStackTrace(void* context = 0);
} // End namespace Uintah

#endif
