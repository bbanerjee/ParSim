#ifndef PARELLIP3D_CORE_UTIL_H
#define PARELLIP3D_CORE_UTIL_H

#include <Core/Types/realtypes.h>
#include <string>
#include <sys/time.h>

namespace util {

struct timeval timediff(const struct timeval& time1,
                        const struct timeval& time2);

long int timediffmsec(const struct timeval& time1, const struct timeval& time2);

REAL timediffsec(const struct timeval& time1, const struct timeval& time2);

std::string combine(const std::string& folder, const std::string& str, std::size_t num, std::size_t width);

template <typename T>
T getParam(const std::string str);

// Creates output folder and returns name
std::string createOutputFolder(const std::string& folderName);

} // end namespace util
#endif // PARELLIP3D_CORE_UTIL_H