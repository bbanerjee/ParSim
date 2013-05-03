#ifndef EMU2DC_EXCEPTION_H
#define EMU2DC_EXCEPTION_H

#include <stdexcept>
#include <string>
#include <sstream>

namespace Emu2DC {

  class Exception : public std::runtime_error
  {
  public:
    
    Exception(const std::string& msg, const char* file, int line):
	std::runtime_error("")
    {
      std::ostringstream s;
      s << "Exception thrown: " << file << ", line: " << line << "\n" << msg;
      static_cast<std::runtime_error&>(*this) = std::runtime_error(s.str());
    }

    Exception(const std::ostringstream& msg, const char* file, int line):
	std::runtime_error("")
    {
      std::ostringstream s;
      s << "Exception thrown: " << file << ", line: " << line << "\n" << msg;
      static_cast<std::runtime_error&>(*this) = std::runtime_error(s.str());
    }
  };

} // end namespace
#endif
