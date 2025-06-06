/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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
 *  Assert.h: Utility for specifying data invariants (Assertions)
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1994
 *
 */

#ifndef SCI_Containers_FancyAssert_h
#define SCI_Containers_FancyAssert_h

#include <sci_defs/error_defs.h>

#include <Core/Exceptions/AssertionFailed.h>
#include <sstream>

/*
 * Note - a normal AssertionFailed exception cannot be used here.  We
 * must use one that takes a string
 */

#if SCI_ASSERTION_LEVEL >= 2
#define ASSERTEQ(c1, c2) \
   if(c1 != c2){ \
      std::ostringstream msg; \
      msg << #c1 << "(value=" << c1 << ") == " << #c2 << "(value=" << c2 << ")"; \
      SCI_THROW(Uintah::AssertionFailed(msg.str().c_str(), __FILE__, __LINE__)); \
   }

#define ASSERTNE(c1, c2) \
   if(c1 != c2){ \
      std::ostringstream msg; \
      msg << #c1 << "(value=" << c1 << ") != " << #c2 << "(value=" << c2 << ")"; \
      SCI_THROW(Uintah::AssertionFailed(msg.str().c_str(), __FILE__, __LINE__)); \
   }

#define ASSERTRANGE(c, l, h) \
   if(c < l || c >= h){ \
      std::ostringstream msg; \
      msg << #l "(value=" << l << ") <= " #c << "(value=" << c << ") < " << #h << "(value=" << h << ")"; \
      SCI_THROW(Uintah::AssertionFailed(msg.str().c_str(), __FILE__, __LINE__)); \
   }
#else
#define ASSERTEQ(c1, c2)
#define ASSERTRANGE(c, l, h)
#endif

#endif
