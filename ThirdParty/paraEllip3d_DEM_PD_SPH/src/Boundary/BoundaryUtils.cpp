#include <Boundary/BoundaryUtils.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>

#include <vector>
#include <string>

namespace dem {

  namespace BCUtils {

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
    Displacement convert<Displacement>(const std::vector<double>& vec)
    {
      return Displacement(vec[0], vec[1], vec[2]);
    }

    template<> 
    DeformationGradient convert<DeformationGradient>(const std::vector<double>& vec)
    {
      return DeformationGradient(vec[0], vec[1], vec[2],
                                 vec[3], vec[4], vec[5],
                                 vec[6], vec[7], vec[8]);
    }

    template<> 
    AxisymmetricStrain convert<AxisymmetricStrain>(const std::vector<double>& vec)
    {
      return AxisymmetricStrain(vec[0], vec[1], vec[2], vec[3]);
    }

    template<> 
    std::vector<double> 
    convertStrToArrayOfScalars<double>(const std::string& vecStr)
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

    /*
    template<>
    std::vector<Vec> convertStrToArray<Vec>(const std::string& vecStr)
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
    */

    template <>
    double 
    createValue<double>(const std::vector<double>& vec, std::size_t index)
    {
      return vec[index];
    }

    template <>
    Vec 
    createValue<Vec>(const std::vector<double>& vec, std::size_t index)
    {
      return Vec(vec[3*index], vec[3*index+1], vec[3*index+2]);
    }

    template <>
    Displacement 
    createValue<Displacement>(const std::vector<double>& vec, std::size_t index)
    {
      return Displacement(vec[3*index], vec[3*index+1], vec[3*index+2]);
    }

    template <>
    DeformationGradient 
    createValue<DeformationGradient>(const std::vector<double>& vec, std::size_t index)
    {
      auto ii = 9*index;
      return DeformationGradient(vec[ii], vec[ii+1], vec[ii+2],
                                 vec[ii+3], vec[ii+4], vec[ii+5],
                                 vec[ii+6], vec[ii+7], vec[ii+8]);
    }

    template <>
    AxisymmetricStrain 
    createValue<AxisymmetricStrain>(const std::vector<double>& vec, std::size_t index)
    {
      auto ii = 4*index;
      return AxisymmetricStrain(vec[ii], vec[ii+1], vec[ii+2], vec[ii+3]);
    }

  } // end namespace BCUtils
}// end namespace dem