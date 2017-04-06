#include <Core/Util/Utility.h>
#include <InputOutput/Parameter.h>
#include <string>
#include <sstream>
#include <iomanip>

namespace util {

struct timeval
timediff(const struct timeval& time1, const struct timeval& time2)
{
  struct timeval diff;
  diff.tv_sec = time2.tv_sec - time1.tv_sec;
  if ((diff.tv_usec = time2.tv_usec - time1.tv_usec) < 0) {
    diff.tv_usec += 1000000;
    --diff.tv_sec;
  }
  return (diff);
}

long int
timediffmsec(const struct timeval& time1, const struct timeval& time2)
{
  struct timeval diff = timediff(time1, time2);
  return (diff.tv_sec * 1000000 + diff.tv_usec);
}

REAL
timediffsec(const struct timeval& time1, const struct timeval& time2)
{
  return ((REAL)timediffmsec(time1, time2) / 1.0e+6);
}

std::string
combine(const std::string& str, std::size_t num, std::size_t width)
{
  std::string out(str); 
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  out += ss.str();
  return out;
}

template <typename T> 
T getParam(const std::string str) {
  return static_cast<T>(dem::Parameter::get().param[str]);
}

template <>
std::size_t getParam(const std::string str) {
  return static_cast<std::size_t>(dem::Parameter::get().param[str]);
}

template <>
int getParam(const std::string str) {
  return static_cast<int>(dem::Parameter::get().param[str]);
}

template <>
REAL getParam(const std::string str) {
  return static_cast<REAL>(dem::Parameter::get().param[str]);
}

} // end namespace util