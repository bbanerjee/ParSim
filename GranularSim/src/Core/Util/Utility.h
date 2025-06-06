#ifndef PARELLIP3D_CORE_UTIL_H
#define PARELLIP3D_CORE_UTIL_H

#include <Core/Types/RealTypes.h>
#include <string>
#include <vector>
#include <sys/time.h>

namespace util {

struct timeval timediff(const struct timeval& time1,
                        const struct timeval& time2);

long int timediffmsec(const struct timeval& time1, const struct timeval& time2);

REAL timediffsec(const struct timeval& time1, const struct timeval& time2);

std::string combine(const std::string& folder, const std::string& str, std::size_t num, std::size_t width);
std::string combine(const std::string& str, std::size_t num, std::size_t width);

template <typename T>
T getParam(const std::string str);

std::string getFilename(const std::string& str);

// Creates output folder and returns name
std::string createOutputFolder(const std::string& folderName);

// Create a std::vector of equally spaced ints/reals
template <typename T>
std::vector<T> linspaceApprox(const T& low, const T& high, T spacing);
template <typename T>
std::vector<T> linspace(const T& low, const T& high, int numSteps);

// Replace an extension (e.g., ".yyy") with (e.g., ."zzz") in a string 
// of the form "xxxx.yyy"
std::string replaceExtension(const std::string& str, const std::string& newExtension);

} // end namespace util
#endif // PARELLIP3D_CORE_UTIL_H