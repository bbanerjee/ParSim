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

namespace Ellip3D {
namespace IOUtil {

template <typename T>
bool
decodeAndUncompress(const std::string& inputStr,
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
    output.push_back(convert<T>(str));
  }

  return true;
}

template <typename T>
bool
compressAndEncode(const std::vector<T>& inputVec,
                                 const int& numComponents,
                                 std::string& outputStr)
{
  // Convert the vector into a string
  std::string inputStr;
  for (const auto& val : inputVec) {
    inputStr += " ";
    inputStr += convert<T>(val);
  }

  // Compress the string
  // See: https://panthema.net/2007/0328-ZLibString.html
  z_stream stream;
  stream.zalloc = Z_NULL;
  stream.zfree = Z_NULL;
  stream.opaque = Z_NULL;

  int err = deflateInit(&stream, Z_BEST_COMPRESSION);
  if (err != Z_OK) {
    std::cerr << "zlib::deflateInit" << " error: " << err 
              << ":" << stream.msg << std::endl;
    return false;
  }

  stream.avail_in = inputStr.size();
  stream.next_in = (Bytef *) inputStr.data();

  char buffer[32768];
  std::string compressedStr;

  do {
    stream.next_out = reinterpret_cast<Bytef *>(buffer);
    stream.avail_out = sizeof(buffer);

    err = deflate(&stream, Z_FINISH);
    if (compressedStr.size() < stream.total_out) {
      compressedStr.append(buffer, stream.total_out - compressedStr.size());

    }
  } while (err == Z_OK);

  deflateEnd(&stream);
  if (err != Z_STREAM_END) {
    std::cerr << "zlib::deflate" << " error: " << err 
              << ":" << stream.msg << std::endl;
    return false;
  }

  // Encode the compressed string
  outputStr = base64::encode(compressedStr);

  return true;
}

template <>
bool
convert<bool>(const std::string& str)
{
  bool val;
  std::istringstream(str) >> std::boolalpha >> val;
  return val;
}

template <>
int
convert<int>(const std::string& str)
{
  return std::stoi(str);
}

template <>
size_t
convert<size_t>(const std::string& str)
{
  return std::stoul(str);
}

template <>
REAL
convert<REAL>(const std::string& str)
{
  return std::stod(str);
}

template <>
Vec
convert<Vec>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return Vec(std::stod(split[0]), std::stod(split[1]), std::stod(split[2]));
}

template <>
IntVec
convert<IntVec>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return IntVec(std::stoi(split[0]), std::stoi(split[1]), std::stoi(split[2]));
}

template <>
std::vector<REAL>
convert<std::vector<REAL>>(const std::string& str)
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  std::vector<REAL> vec;
  for (auto str : split) {
    vec.push_back(std::stod(str));
  }
  return vec;
}

template <>
std::string 
convert<std::size_t>(const std::size_t& value)
{
  return std::to_string(value);
}

template <>
std::string 
convert<int>(const int& value)
{
  return std::to_string(value);
}

template <>
std::string 
convert<REAL>(const REAL& value)
{
  return std::to_string(value);
}

template <>
std::string 
convert<Vec>(const Vec& value)
{
  return std::to_string(value.x()) + " " +
         std::to_string(value.y()) + " " + 
         std::to_string(value.z());
}

template <>
std::string 
convert<IntVec>(const IntVec& value)
{
  return std::to_string(value.x()) + " " +
         std::to_string(value.y()) + " " + 
         std::to_string(value.z());
}

template <>
std::string 
convert<std::vector<REAL>>(const std::vector<REAL>& value)
{
  std::string out;
  for (const auto& val : value) {
    out += std::to_string(val) + " ";
  }
  return out;
}

template <typename T>
std::vector<T>
convertStrArray(const std::string& str)
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
decodeAndUncompress<int>(const std::string& inputStr,
  const int& numComponents, std::vector<int>& output);

template bool
decodeAndUncompress<std::size_t>(const std::string& inputStr,
  const int& numComponents, std::vector<std::size_t>& output);

template bool
decodeAndUncompress<REAL>(const std::string& inputStr,
  const int& numComponents, std::vector<REAL>& output);

template bool
decodeAndUncompress<dem::IntVec>(const std::string& inputStr,
  const int& numComponents, std::vector<dem::IntVec>& output);

template bool
decodeAndUncompress<dem::Vec>(const std::string& inputStr,
  const int& numComponents, std::vector<dem::Vec>& output);

template bool
decodeAndUncompress<std::vector<REAL>>(const std::string& inputStr,
  const int& numComponents, std::vector<std::vector<REAL>>& output);

template bool
compressAndEncode<std::size_t>(const std::vector<std::size_t>& inputVec,
  const int& numComponents, std::string& outputStr);

template bool
compressAndEncode<int>(const std::vector<int>& inputVec,
  const int& numComponents, std::string& outputStr);

template bool
compressAndEncode<REAL>(const std::vector<REAL>& inputVec,
  const int& numComponents, std::string& outputStr);

template bool
compressAndEncode<IntVec>(const std::vector<IntVec>& inputVec,
  const int& numComponents, std::string& outputStr);

template bool
compressAndEncode<Vec>(const std::vector<Vec>& inputVec,
  const int& numComponents, std::string& outputStr);

template bool
compressAndEncode<std::vector<REAL>>(const std::vector<std::vector<REAL>>& inputVec,
  const int& numComponents, std::string& outputStr);

template 
std::vector<REAL>
convertStrArray<REAL>(const std::string& str);

} // end namespace IOUtil
} // end namespace Ellip3D
