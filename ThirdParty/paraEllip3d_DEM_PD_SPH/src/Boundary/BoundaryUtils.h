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
    S createValue(const std::vector<double>& vec, std::size_t index);
    
    template <typename S> 
    std::vector<S> convertStrToArrayOfScalars(const std::string& vecStr);

    template <typename S, int N> 
    std::vector<S> convertStrToArrayOfVectors(const std::string& vecStr)
    {
      std::istringstream iss(std::string(vecStr.begin(), vecStr.end()));
      std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                          std::istream_iterator<std::string>{} };
      std::vector<double> vec;
      for (auto str : split) {
        vec.push_back(std::stod(str));
      }
      std::vector<S> valVec;
      valVec.resize(vec.size()/N);
      for (auto ii = 0u; ii < valVec.size(); ii++) {
        //valVec[ii] = S(vec[3*ii], vec[3*ii+1], vec[3*ii+2]);
        valVec[ii] = createValue<S>(vec, ii);
      }
      return valVec;
    }

    /*
    template<> 
    double convert<double>(const std::vector<double>& vec);

    template<> 
    Vec convert<Vec>(const std::vector<double>& vec);

    template<> 
    std::vector<double> convertStrToArrayOfScalars<double>(const std::string& vecStr);

    template<>
    std::vector<Vec> convertStrToArrayOfVectors<Vec, 3>(const std::string& vecStr);
    */

  } // end namespace BCUtils
}
#endif