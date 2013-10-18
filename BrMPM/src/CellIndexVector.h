/*
 * CellIndexVector.h
 *
 *  Created on: 18/10/2013
 *      Author: banerjee
 */

#ifndef CELLINDEXVECTOR_H_
#define CELLINDEXVECTOR_H_

#include <vector>

namespace BrMPM {

  class CellIndexVector
  {
  public:

    CellIndexVector() {}
    virtual ~CellIndexVector() {}

    inline CellIndexVector(const CellIndexVector& vec)
    {
      for (auto iter = (vec.d_data).begin(); iter != (vec.d_data).end(); ++iter) {
        d_data.emplace_back(*iter);
      }
    }

    inline CellIndexVector(const std::vector<int>& data)
    {
      for (auto iter = data.begin(); iter != data.end(); ++iter) {
        d_data.emplace_back(*iter);
      }
    }

    inline void operator=(const CellIndexVector& vec) {
      for (unsigned int ii = 0; ii < (vec.d_data).size(); ++ii) {
        d_data[ii] = (vec.d_data)[ii];
      }
    }

    inline unsigned int size() const {return d_data.size();}

    inline int& operator[](int index) {return d_data.at(index);}
    inline const int& operator[](int index) const {return d_data.at(index);}

    inline CellIndexVector operator*(const CellIndexVector& vec) const {
      return CellIndexVector(*this);
    }

  private:

    std::vector<int> d_data;
  };

} /* namespace BrMPM */

#endif /* CELLINDEXVECTOR_H_ */
