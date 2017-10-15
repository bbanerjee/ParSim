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
    std::vector<S> convertStrToArray(const std::string& vecStr);
    
    template<> 
    double convert<double>(const std::vector<double>& vec);

    template<> 
    Vec convert<Vec>(const std::vector<double>& vec);

    template<> 
    std::vector<double> convertStrToArray<double>(const std::string& vecStr);

    template<>
    std::vector<Vec> convertStrToArray<Vec>(const std::string& vecStr);

  } // end namespace BCUtils
}
#endif