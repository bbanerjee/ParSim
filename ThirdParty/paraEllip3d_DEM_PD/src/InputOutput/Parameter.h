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

class Parameter
{

public:
  // static function is part of the class, not part of the object, so it is used
  // like Parameter::getSingleton(). it can only access static members.
  // it is public so that it can be called by others like
  // dem::Parameter::getSingleton()
  // and it can implicitly call private constructor Parameter().
  static Parameter& getSingleton()
  {
    static Parameter instance; // instantiated on first use,
                               // guaranteed to be destroyed
    return instance;
  }

  void readIn(const char* input);
  int  readInXML(const std::string& inputFileName);
  void writeOut();
  void writeOutXML();

private:
  // constructor must be private to avoid instantiation by others because
  // singleton
  // should only be instantiated by itself
  Parameter() {}
  ~Parameter() {}
  // make sure these two are unaccessable to avoid copies of singelton
  Parameter(Parameter const&) = delete;      // don't implement
  void operator=(Parameter const&) = delete; // don't implement

public:
  std::map<std::string, REAL> parameter;
  std::vector<std::pair<REAL, REAL>> gradation;
  std::map<std::string, std::string> datafile;
  std::vector<REAL> sigmaPath;

private:
  friend class boost::serialization::access;
  template <class ArchiveType>
  void serialize(ArchiveType& ar, const unsigned int version)
  {
    ar& parameter;
    ar& sigmaPath;
  }
};
}
#endif
