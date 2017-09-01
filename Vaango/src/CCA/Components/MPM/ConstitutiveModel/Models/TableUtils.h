#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H

#include <string>
#include <vector>

namespace Vaango {

namespace Util {

// Taken from:
// https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string/236803#236803
template <typename Out>
void split(const std::string& s, char delim, Out result);
std::vector<std::string> split(const std::string& s, char delim);

// Taken from:
// https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v");
std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v");
std::string& trim(std::string& s, const char* t = " \t\n\r\f\v");

double toDouble(std::string& str);
}
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H
