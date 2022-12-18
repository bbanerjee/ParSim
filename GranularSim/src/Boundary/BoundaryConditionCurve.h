#ifndef ELLIP3D_BOUNDARY_CONDITION_CURVE_H
#define ELLIP3D_BOUNDARY_CONDITION_CURVE_H

#include <Boundary/BoundaryUtils.h>
#include <InputOutput/json/json.hpp>
#include <InputOutput/zenxml/xml.h>
#include <Core/Math/Vec.h>
#include <Core/MechanicsConcepts/Deformations.h>

#include <vector>
#include <fstream>
#include <algorithm>

namespace dem {

  using XMLProblemSpec = zen::XmlIn;
  using JsonProblemSpec = nlohmann::json;

  template<typename T, int N>
  class BoundaryConditionCurve {

  public:
    BoundaryConditionCurve() = default;
    BoundaryConditionCurve(const BoundaryConditionCurve&) = delete;
    BoundaryConditionCurve& operator=(const BoundaryConditionCurve&) = delete;
    ~BoundaryConditionCurve() = default;

    BoundaryConditionCurve(std::istream& csv_format)
    {
      read(csv_format);
    }

    BoundaryConditionCurve(const XMLProblemSpec& xml_format)
    {
      read(xml_format);
    }

    BoundaryConditionCurve(const JsonProblemSpec& json_format)
    {
      read(json_format);
    }

    void
    read(std::istream& ifs)
    {
      int num_points = 0, num_components = 0;
      ifs >> num_points >> num_components;

      d_times.resize(num_points);
      d_bcValues.resize(num_points);

      for (int ii = 0; ii < num_points; ii++) {
        ifs >> d_times[ii];
      }
      for (int ii = 0; ii < num_points; ii++) {
        std::vector<double> valueVec(num_components);
        for (int jj = 0; jj < num_components; jj++) {
          ifs >> valueVec[jj];
        }
        d_bcValues[ii] = dem::BCUtils::convert<T>(valueVec);
      }
      assert(d_times.size() == d_bcValues.size());
    }

    void
    read(const XMLProblemSpec& ps)
    {
      std::string vecStr;
      if (!ps["time"](vecStr)) {
        std::cerr
          << "**ERROR** Time data not found in boundary condition curve\n";
        std::cerr << "  Add the <time> t1, t2 t3, ... </time> tag.";
        exit(-1);
      }
      d_times = dem::BCUtils::convertStrToArrayOfScalars<double>(vecStr);

      vecStr = "";
      if (!ps["value"](vecStr)) {
        std::cerr
          << "**ERROR** Load/displacement data not found in boundary condition curve\n";
        std::cerr << "  Add the <value> v1, v2 v3, ... </value> tag.";
        exit(-1);
      }
      d_bcValues = dem::BCUtils::convertStrToArrayOfVectors<T, N>(vecStr);
      assert(d_times.size() == d_bcValues.size());
    }

    void
    read(const JsonProblemSpec& ps)
    {
      std::string vecStr;
      try {
        vecStr = ps["time"].get<std::string>();
      } catch (std::exception& e) {
        std::cerr
          << "**ERROR** Time data not found in boundary condition curve\n";
        std::cerr << "  Add the time: \"t1 t2 t3 ...\" tag.";
        exit(-1);
      }
      d_times = dem::BCUtils::convertStrToArrayOfScalars<double>(vecStr);
    
      vecStr = "";
      try {
        vecStr = ps["value"].get<std::string>();
      } catch (std::exception& e) {
        std::cerr
          << "**ERROR** Load/displacement data not found in boundary condition curve\n";
        std::cerr << "  Add the value: \"v1 v2 v3 ...\" tag.";
        exit(-1);
      }
      d_bcValues = dem::BCUtils::convertStrToArrayOfVectors<T, N>(vecStr);
      assert(d_times.size() == d_bcValues.size());
    }


    /**
     * Linear interpolation
     * REQUIRES: time is monotonically increasing
     */
    T getBCValue(double time) const {

      if (d_times.size() < 2 || time <= d_times.front()) return d_bcValues.front();
      if (time >= d_times.back()) return d_bcValues.back();

      auto iter_hi = std::upper_bound(d_times.begin(), d_times.end(), time);
      auto index = iter_hi - d_times.begin();
      double s = (time  - d_times[index-1])/(d_times[index] - d_times[index-1]);
      return (1 - s)*d_bcValues[index-1] + s*d_bcValues[index];
    }

    /**
     * Time Derivative
     * REQUIRES: time is monotonically increasing
     */
    T getBCRate(double time) const {

      if (d_times.size() < 2 || time < d_times.front() || time > d_times.back()) {
        return T(0);
      }

      double lo = time - 1.0e-6;
      double hi = time + 1.0e-6;
      if (lo < d_times.front()) lo = d_times.front();
      if (hi > d_times.back()) hi = d_times.back();

      auto lo_val = getBCValue(lo);
      auto hi_val = getBCValue(hi);
      return (hi_val - lo_val)*2.0e6;
    }

    friend std::ostream& operator<<(std::ostream& os, const BoundaryConditionCurve& bc)
    {
      os << "Times: ";
      for (auto time : bc.d_times) {
        os << time << " ";
      }
      os << std::endl;
      os << "Values: ";
      for (auto value : bc.d_bcValues) {
        os << value << " ";
      }
      os << std::endl;
      return os;
    }

  private:
    std::vector<double> d_times;
    std::vector<T> d_bcValues;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      ar& d_times;
      ar& d_bcValues;
    }

  };

} // end namespace dem

namespace dem {
  template class BoundaryConditionCurve<double, 1>;
  template class BoundaryConditionCurve<Vec, 3>;
  template class BoundaryConditionCurve<Displacement, 3>;
}

#endif