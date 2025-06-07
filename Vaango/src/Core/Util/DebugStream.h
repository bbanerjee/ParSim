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

/*! DebugStream.h - An ostream used for debug messages
 *
 *  Written by:
 *  Eric Kuehne
 * Department of Computer Science
 * University of Utah
 * Feb. 2000
 *
 * DebugStream is an ostream that is useful for outputing debug messages.
 * When an instance is created, it is given a name.  An environment variable,
 * SCI_DEBUG, is inspected to see if a particular instance should be
 * active(identified by its name), and if so where to send the output.
 * The syntax for the environment variable is:
 * SCI_DEBUG = ([name]:[-|+|+FILENAME])(,[name]:[-|+|+FILENAME])*
 * The + or - specifies wheather the named object is on or off.  If a file is
 * specified it is opened in ios::out mode.  If no file is specified,
 * the stream is directed to cerr.  The : and , characters are
 * restricted to deliminators.
 *
 * Example:
 * SCI_DEBUG = modules.meshgen.warning:+meshgen.out,util.debugstream.error:-
 *
 * Future Additions:
 * o Possible additions to constructor:
 *   - Default file to output to
 *   - Mode that the file will be opened in (append, out, etc.) (usefulness??)
 * o Allow DEFAULT specification in the env variable which would override
 *   all default settings. (usefulness??)
 * o Time stamp option
 *
 * Annoyances:
 * o Because the list of environment variables, "environ", is built at
 * run time, and getenv queries that list, I have not been able to
 * figure out a way to requery the environment variables during
 * execution.
 */

#ifndef __VAANGO_CORE_UTIL_DebugStream_h__
#define __VAANGO_CORE_UTIL_DebugStream_h__

#include <iomanip>
#include <iostream>
#include <map>
#include <string>

namespace Uintah {

class DebugStream;
class DebugBuf;

///////////////////
// class DebugBuf
// For use with DebugStream.  This class overrides the overflow
// operator.  Each time overflow is called it checks to see where
// to direct the output to.
class DebugBuf : public std::streambuf
{
private:
public:
  DebugBuf() = default;
  ~DebugBuf() = default;

  int
  overflow(int ch) override;

  // points the the DebugStream that instantiated me
  DebugStream* owner{ nullptr };
};

///////////////////
// class DebugStream
// A general purpose debugging ostream.
class DebugStream : public std::ostream
{
public:
  inline static std::map<std::string, DebugStream*> s_all_debug_streams{};
  inline static bool s_all_dbg_streams_initialized{ false };

public:
  DebugStream();

  DebugStream(const std::string& name, bool defaulton = true);

  DebugStream(const std::string& name,
              const std::string& component,
              const std::string& description,
              bool defaulton = true);

  ~DebugStream();

  std::string
  getName()
  {
    return d_name;
  }

  std::string
  getComponent()
  {
    return d_component;
  }

  std::string
  getDescription()
  {
    return d_description;
  }

  std::string
  getFilename()
  {
    return d_filename;
  }

  void
  setFilename(const std::string& name)
  {
    d_filename = name;
  }

  bool
  active()
  {
    return d_active;
  };

  void
  setActive(bool active)
  {
    d_active = active;
  };

  void
  print() const
  {
    std::cout << std::setw(2) << std::left << (d_active ? "+" : "-")
              << std::setw(30) << std::left << d_name.c_str() << std::setw(75)
              << std::left << d_description.c_str() << std::setw(30)
              << std::left << d_component.c_str() << std::endl;
  }

  static void
  printAll()
  {
    printf("-------------------------------------------------------------------"
           "-------------\n");
    for (auto& [name, stream] : s_all_debug_streams) {
      stream->print();
    }
    printf("-------------------------------------------------------------------"
           "-------------\n\n");
  }

  // the ostream that output should be redirected to. std::cout by default.
  std::ostream* outstream{ nullptr };

private:
  // identifies me uniquely
  std::string d_name{ "" };
  std::string d_component{ "Unknown" };
  std::string d_description{ "No description" };

  // if false, all input is ignored
  bool d_active{ false };

  // The stream filename.
  std::string d_filename{ "cout" };

  // the buffer that is used for output redirection
  DebugBuf d_dbgbuf{};

  // check the environment variable
  void
  checkEnv();

  // Check for a previous name.
  void
  checkName();

  void
  instantiate_map();
};

} // End namespace Uintah

#endif //__VAANGO_CORE_UTIL_DebugStream_h__
