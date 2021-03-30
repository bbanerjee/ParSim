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

#ifndef MATITI_EXCEPTION_H
#define MATITI_EXCEPTION_H

#include <stdexcept>
#include <string>
#include <sstream>

namespace Matiti {

  class Exception : public std::runtime_error
  {
  public:
    
    Exception(const std::string& msg, const char* file, int line):
	std::runtime_error("")
    {
      std::ostringstream s;
      s << "Exception thrown: " << file << ", line: " << std::to_string(line) << "\n" << msg;
      static_cast<std::runtime_error&>(*this) = std::runtime_error(s.str());
    }

    Exception(const std::ostringstream& msg, const char* file, int line):
	std::runtime_error("")
    {
      std::ostringstream s;
      s << "Exception thrown: ";
      s << file << ", line: ";
      s << std::to_string(line);
      s << "\n";
      s << msg.str();
      static_cast<std::runtime_error&>(*this) = std::runtime_error(s.str());
    }
  };

} // end namespace
#endif