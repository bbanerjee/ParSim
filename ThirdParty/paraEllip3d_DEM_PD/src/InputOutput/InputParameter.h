/*
In general, built-in C++ types (ints, floats, characters, etc.) can be
transmitted over MPI directly.

For types defined by the standard library (such as std::string or std::vector)
and some types in Boost
(such as boost::variant), the Boost.Serialization library already contains all
of the required
serialization code. In these cases, you need only include the appropriate header
from the boost/serialization
directory.

For types that do not already have a serialization header, you will first need
to implement serialization
code before the types can be transmitted using Boost.MPI.
*/

#ifndef PARAMETER_H
#define PARAMETER_H

#include <Core/Types/realtypes.h>
#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace dem {

class InputParameter
{

public:
  // static function is part of the class, not part of the object, so it is used
  // like InputParameter::get(). it can only access static members.
  // it is public so that it can be called by others like
  // InputParameter::get()
  // and it can implicitly call private constructor InputParameter().
  static InputParameter& get()
  {
    static InputParameter instance; // instantiated on first use,
                               // guaranteed to be destroyed
    return instance;
  }

  void readIn(const char* input);
  bool readInXML(const std::string& inputFileName);
  void writeOut();
  void writeOutXML();
  void addParameter(const std::string& key, REAL value) {
    param[key] = value;
  }

private:
  // constructor must be private to avoid instantiation by others because
  // singleton
  // should only be instantiated by itself
  InputParameter() = default;
  ~InputParameter() = default;
  // make sure these two are unaccessable to avoid copies of singelton
  InputParameter(InputParameter const&) = delete;      // don't implement
  void operator=(InputParameter const&) = delete; // don't implement

public:
  std::map<std::string, REAL> param;
  std::vector<std::pair<REAL, REAL>> gradation;
  std::map<std::string, std::string> datafile;
  std::vector<REAL> sigmaPath;

private:
  friend class boost::serialization::access;
  template <class ArchiveType>
  void serialize(ArchiveType& ar, const unsigned int version)
  {
    ar& param;
    ar& sigmaPath;
  }
};
}
#endif
