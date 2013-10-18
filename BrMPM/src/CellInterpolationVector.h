/*
 * CellInterpolationVector.h
 *
 *  Created on: 18/10/2013
 *      Author: banerjee
 */

#ifndef CELLINTERPOLATIONVECTOR_H_
#define CELLINTERPOLATIONVECTOR_H_

namespace BrMPM {

  class CellInterpolationVector
  {
  public:
    CellInterpolationVector() {}
    virtual ~CellInterpolationVector() {}

    CellInterpolationVector(const CellInterpolationVector& vec)
    {
      for (auto iter = (vec.d_data).begin(); iter != (vec.d_data).end(); ++iter) {
        d_data.emplace_back(*iter);
      }
    }

    CellInterpolationVector(const std::vector<double>& data)
    {
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        d_data.emplace_back(*iter);
      }
    }

    void operator=(const CellInterpolationVector& vec) {
      for (unsigned int ii = 0; ii < (vec.d_data).size(); ++ii) {
        d_data[ii] = (vec.d_data)[ii];
      }
    }

    unsigned int size() const {return d_data.size();}

    inline double& operator[](int index) {return d_data.at(index);}
    inline const double& operator[](int index) const {return d_data.at(index);}

    inline CellInterpolationVector operator*(const CellInterpolationVector& vec) const {
      return CellInterpolationVector(*this);
    }

  private:

    std::vector<double> d_data;
  };

} /* namespace BrMPM */
#endif /* CELLINTERPOLATIONVECTOR_H_ */
