#ifndef ELLIP3D_BOUNDARY_UTILS_H
#define ELLIP3D_BOUNDARY_UTILS_H

#include <Core/Math/Vec.h>

#include <vector>
#include <string>

namespace dem {

  namespace BCUtils {

    template <typename S> 
    S convert(const std::vector<double>& vec);

    template <typename S> 
    std::vector<S> convert(const std::string& vecStr);
    
    template<> 
    double convert<double>(const std::vector<double>& vec)
    {
      return vec[0];
    }

    template<> 
    Vec convert<Vec>(const std::vector<double>& vec)
    {
      return Vec(vec[0], vec[1], vec[2]);
    }

    template<> 
    std::vector<double> convert<double>(const std::string& vecStr)
    {
      std::istringstream iss(std::string(vecStr.begin(), vecStr.end()));
      std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                          std::istream_iterator<std::string>{} };
      std::vector<double> vec;
      for (auto str : split) {
        vec.push_back(std::stod(str));
      }
      return vec;
    }

    template<> 
    std::vector<Vec> convert<Vec>(const std::string& vecStr)
    {
      std::istringstream iss(std::string(vecStr.begin(), vecStr.end()));
      std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                          std::istream_iterator<std::string>{} };
      std::vector<double> vec;
      for (auto str : split) {
        vec.push_back(std::stod(str));
      }
      std::vector<Vec> valVec;
      valVec.resize(vec.size()/3);
      for (auto ii = 0u; ii < valVec.size(); ii++) {
        valVec[ii] = Vec(vec[3*ii], vec[3*ii+1], vec[3*ii+2]);
      }
      return valVec;
    }
  } // end namespace BCUtils
}
#endif