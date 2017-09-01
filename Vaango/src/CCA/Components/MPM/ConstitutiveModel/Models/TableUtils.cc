#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>
#include <algorithm>
#include <sstream>

namespace Vaango {

namespace Util {

// Taken from:
// https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string/236803#236803
template <typename Out>
void
split(const std::string& s, char delim, Out result)
{
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string>
split(const std::string& s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

// trim from left
std::string&
ltrim(std::string& s, const char* t)
{
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from right
std::string&
rtrim(std::string& s, const char* t)
{
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from left & right
std::string&
trim(std::string& s, const char* t)
{
  return ltrim(rtrim(s, t), t);
}

// Convert into double
double
toDouble(std::string& str)
{
  // Remove illegal characters
  const std::string illegalChars = " [\\/:?\"<>|]";
  str.erase(remove_if(str.begin(), str.end(),
                      [&illegalChars](char cc) {
                        return (illegalChars.find(cc) != std::string::npos)
                                 ? true
                                 : false;
                      }),
            str.end());

  // Convert the string to double
  return std::stod(str);
}

} // End namespace Util

} // end namespace Vaango