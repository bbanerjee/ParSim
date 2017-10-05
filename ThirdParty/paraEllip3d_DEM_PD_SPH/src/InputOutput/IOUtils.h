#ifndef ELLIP3D_INPUTOUTPUT_IO_UTIL_H
#define ELLIP3D_INPUTOUTPUT_IO_UTIL_H

#include <vector>
#include <string>

namespace Ellip3D {

namespace Util {

  /* A reverse range iterator 
    std::vector<int> v = {1, 2, 3, 4, 5};
    for (auto x : reverse(v)) {}
    From : https://gist.github.com/arvidsson/7231973 */
  template <typename T>
  class reverse_range
  {
    T &x;

  public:
    reverse_range(T &x) : x(x) {}
    auto begin() const -> decltype(this->x.rbegin())
    {
      return x.rbegin();
    }
    auto end() const -> decltype(this->x.rend())
    {
      return x.rend();
    }
  };
  template <typename T>
  reverse_range<T> reverse(T &x)
  {
    return reverse_range<T>(x);
  }

  /**
   * Decode and uncompress a base64 string
   */
  template <typename T>
  bool decodeAndUncompress(const std::string& inputStr,
                           const int& numComponents,
                           std::vector<T>& output);

  /**
   * Compress a string created from a std::vector and then encode to base64
   */
  template <typename T>
  bool compressAndEncode(const std::vector<T>& input,
                         const int& numComponents,
                         std::string& outputStr);

  /**
   * Convert a string to another type
   */
  template <typename T>
  T convert(const std::string& str);

  /**
   * Convert a value to a string
   */
  template <typename T>
  std::string convert(const T& value);

  /**
   * Convert a string to an array of another type
   */
  template <typename T>
  std::vector<T> convertStrArray(const std::string& str);

} // end namespace Util
} // end namespace Ellip3D

#endif // ELLIP3D_INPUTOUTPUT_IO_UTIL_H
