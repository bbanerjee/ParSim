#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <iomanip>
#include <sstream>
#include <string>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

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
combine(const std::string& folder, const std::string& str, std::size_t num, std::size_t width)
{
  std::string out(folder);
  out = out + "/" + str;
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  out += ss.str();
  return out;
}

std::string 
createOutputFolder(const std::string& folderName)
{
  // Get the path
  std::string currentWorkingDir(".");
  char buffer[2000];
  char* str = getcwd(buffer, 2000);
  if (str == nullptr) {
    std::cerr << "**ERROR** Current working directory not returned by getcwd()" << __FILE__
              << __LINE__ << "\n";
  } else {
    currentWorkingDir = std::string(buffer);
  }
  int dircount = 0;
  std::ostringstream folderNameStream;
  folderNameStream << currentWorkingDir << "/";
  folderNameStream << folderName;
  folderNameStream << "." << std::setfill('0') << std::setw(3) << dircount;

#if defined(_WIN32)
    _mkdir((folderNameStream.str()).c_str());
#else
    while (opendir((folderNameStream.str()).c_str())) {
      ++dircount;
      folderNameStream.clear();
      folderNameStream.str("");
      folderNameStream << folderName;
      folderNameStream << "." << std::setfill('0') << std::setw(3) << dircount;
    }
    mkdir((folderNameStream.str()).c_str(), 0777);
#endif

  return folderNameStream.str();
}

template <typename T>
T
getParam(const std::string str)
{
  auto param = dem::InputParameter::get().param;
  REAL val;
  try {
    val = param.at(str);
  } catch (const std::out_of_range& err) {
    std::cerr << "Required InputParameter [" << str << "] not found in input file" << '\n';
  }
  return static_cast<T>(param[str]);
}

// Create a std::vector of equally spaced ints/reals; if equality with high
// cannot be achieved add the high value to the end
template <typename T>
std::vector<T> linspaceApprox(const T& low, const T& high, T spacing)
{
  int numSpaces = (spacing <= 0) ? 1 : round((high - low)/spacing);
  int numPoints = (numSpaces <= 1) ? 2 : numSpaces+1;
  std::vector<T> output(numPoints, static_cast<T>(0));
  if (numPoints == 2) {
    output[0] = low;
    output[1] = high;
  } else {
    for (int ii = 0; ii < numPoints; ii++) {
      output[ii] = low + static_cast<T>(ii*spacing);
    }
    /*
    if (output[numPoints-1] < high) {
      output.push_back(high);
    }
    */
  }
  return output;
}

// Create a std::vector of equally spaced ints/reals (no approximations)
template <typename T>
std::vector<T> linspace(const T& low, const T& high, int numSpaces)
{
  numSpaces = (numSpaces <= 1) ? 1 : numSpaces;
  T spacing = (high - low)/static_cast<T>(numSpaces);
  spacing = (spacing == 0) ? (high - low) : spacing;
  int numPoints = (high - low)/spacing + 1;

  std::vector<T> output(numPoints, static_cast<T>(0));
  for (int ii = 0; ii < numPoints; ii++) {
    output[ii] = low + static_cast<T>(ii*spacing);
  }
  return output;
}

} // end namespace util

// Specializations
namespace util {

  template std::size_t getParam<std::size_t>(const std::string str);
  template int getParam<int>(const std::string str);
  template REAL getParam<REAL>(const std::string str);

  template std::vector<int> linspaceApprox(const int& low, const int& high, int spacing);
  template std::vector<REAL> linspaceApprox(const REAL& low, const REAL& high, REAL spacing);
  template std::vector<int> linspace(const int& low, const int& high, int numSpaces);
  template std::vector<REAL> linspace(const REAL& low, const REAL& high, int numSpaces);

} // end namespace util
