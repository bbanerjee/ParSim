#include <InputOutput/IOUtils.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>

#include "zlib.h"
#include <cppcodec/cppcodec/base64_default_rfc4648.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>

using IntVec = dem::IntVec;
using Vec = dem::Vec;

template <typename T>
bool
Ellip3D::Util::decodeAndUncompress(const std::string& inputStr,
                                   const int& numComponents,
                                   std::vector<T>& output)
{
  // Decode from base64
  std::vector<std::uint8_t> decoded = base64::decode(inputStr);

  // Uncompress from gzip
  std::vector<std::uint8_t> uncompressed;
  z_stream stream;

  // Allocate inflate state
  stream.zalloc = Z_NULL;
  stream.zfree = Z_NULL;
  stream.opaque = Z_NULL;
  stream.avail_in = 0;
  stream.next_in = Z_NULL;
  int err = inflateInit(&stream);
  if (err != Z_OK) {
    std::cerr << "inflateInit"
              << " error: " << err << std::endl;
    return false;
  }

  // Uncompress until stream ends
  stream.avail_in = decoded.size();
  stream.next_in = &decoded[0];
  do {
    do {
      std::vector<std::uint8_t> out(decoded.size());
      stream.avail_out = out.size();
      stream.next_out = &out[0];

      err = inflate(&stream, Z_SYNC_FLUSH);
      // auto have = decoded.size() - stream.avail_out;
      uncompressed.insert(std::end(uncompressed), std::begin(out),
                          std::end(out));

    } while (stream.avail_out == 0);
  } while (err != Z_STREAM_END);

  if (inflateEnd(&stream) != Z_OK) {
    std::cerr << "inflateEnd"
              << " error: " << err << std::endl;
    return false;
  }

  // Split the uncompressed string into a vector of tokens
  // (Assume that data are space separated)
  // (See: https://stackoverflow.com/questions/236129/split-a-string-in-c)
  std::istringstream iss(std::string(uncompressed.begin(), uncompressed.end()));
  std::vector<std::string> outputStr = { std::istream_iterator<std::string>{
                                           iss },
                                         std::istream_iterator<std::string>{} };
  for (auto iter = outputStr.begin(); iter != outputStr.end();
       iter += numComponents) {
    std::string str = *iter;
    for (int ii = 1; ii < numComponents; ii++) {
      str += " ";
      str += *(iter + ii);
    }
    output.push_back(Ellip3D::Util::convert<T>(str));
  }

  return true;
}

template <>
int
Ellip3D::Util::convert<int>(const std::string& str)
{
  return std::stoi(str);
}

template <>
size_t
Ellip3D::Util::convert<size_t>(const std::string& str)
{
  return std::stoul(str);
}

template <>
double
Ellip3D::Util::convert<double>(const std::string& str)
{
  return std::stod(str);
}

template <>
Vec
Ellip3D::Util::convert<Vec>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return Vec(std::stod(split[0]), std::stod(split[1]), std::stod(split[2]));
}

template <>
IntVec
Ellip3D::Util::convert<IntVec>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return IntVec(std::stoi(split[0]), std::stoi(split[1]), std::stoi(split[2]));
}

template <>
std::vector<double>
Ellip3D::Util::convert<std::vector<double>>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  std::vector<double> vec;
  for (auto str : split) {
    vec.push_back(std::stod(str));
  }
  return vec;
}

template <typename T>
std::vector<T>
Ellip3D::Util::convertStrArray(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  std::vector<T> vec;

  if (std::is_same<T, REAL>::value) {
    for (auto str : split) {
      vec.push_back(std::stod(str));
    }
  } else if (std::is_same<T, size_t>::value) {
    for (auto str : split) {
      vec.push_back(std::stoul(str));
    }
  } else {
    std::cerr << "**ERROR** Conversion of string array is allowed only for"
              << " numeric types\n";
  }
  return vec;
}

template bool
Ellip3D::Util::decodeAndUncompress<int>(const std::string& inputStr,
  const int& numComponents, std::vector<int>& output);

template bool
Ellip3D::Util::decodeAndUncompress<std::size_t>(const std::string& inputStr,
  const int& numComponents, std::vector<std::size_t>& output);

template bool
Ellip3D::Util::decodeAndUncompress<double>(const std::string& inputStr,
  const int& numComponents, std::vector<double>& output);

template bool
Ellip3D::Util::decodeAndUncompress<dem::IntVec>(const std::string& inputStr,
  const int& numComponents, std::vector<dem::IntVec>& output);

template bool
Ellip3D::Util::decodeAndUncompress<dem::Vec>(const std::string& inputStr,
  const int& numComponents, std::vector<dem::Vec>& output);

template bool
Ellip3D::Util::decodeAndUncompress<std::vector<double>>(const std::string& inputStr,
  const int& numComponents, std::vector<std::vector<double>>& output);

template 
std::vector<double>
Ellip3D::Util::convertStrArray<double>(const std::string& str);