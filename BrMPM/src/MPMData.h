/*
 * MPMData.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMDATA_H_
#define MPMDATA_H_

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <Matrix3D.h>
#include <CellIndexVector.h>
#include <CellInterpolationVector.h>

#include <vector>
#include <string>
#include <memory>

namespace BrMPM {

  template <class T> class MPMData {

    public:

      MPMData() {}

      MPMData(int size, const T& value) {
        d_data.resize(size);
        for (int ii = 0; ii < size; ++ii) {
          d_data[ii] = value;
        }
      }

      MPMData(const MPMData<T>& data) {
        d_data.resize(data.size());
        for (unsigned int ii = 0; ii < data.size(); ++ii) {
          d_data[ii] = data[ii];
        }
      }

      MPMData(const std::vector<T>& data) {
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
          d_data.emplace_back(*iter);
        }
      }

      virtual ~MPMData() {}

      void operator=(const MPMData<T>& data) {
        for (unsigned int ii = 0; ii < data.size(); ++ii) {
          d_data[ii] = data[ii];
        }
      }

      typename std::vector<T>::iterator begin() {return d_data.begin();}
      typename std::vector<T>::iterator end() {return d_data.end();}
      typename std::vector<T>::const_iterator begin() const {return d_data.begin();}
      typename std::vector<T>::const_iterator end() const {return d_data.end();}

      unsigned int size() {return d_data.size();}
      unsigned int size() const {return d_data.size();}

      T& operator[](int index) {return d_data.at(index);}
      const T& operator[](int index) const {return d_data.at(index);}

      inline MPMData<T> operator*(const MPMData<T>& data) const {
        std::vector<T> product_vec;
        for (auto iter = d_data.begin(); iter != d_data.end(); ++iter) {
          unsigned int index = iter - d_data.begin();
          product_vec.push_back(d_data[index]*data[index]);
        }
        return MPMData<T>(product_vec);
      }

      inline void zero() {
        for (auto iter = d_data.begin(); iter != d_data.end(); ++iter) {
          d_data[iter-d_data.begin()] *= T(0.0);
        }
      }

      inline void zeroVector() {
        for (auto iter = d_data.begin(); iter != d_data.end(); ++iter) {
          auto vector = *iter;
          for (auto viter = vector.begin(); viter != vector.end(); ++viter) {
            d_data[iter-d_data.begin()][viter-vector.begin()] *= T(0.0);
          }
        }
      }

      inline MPMData<T> clone() {return MPMData<T>(*this);}

      inline void push_back(const T& val) {
        d_data.emplace_back(val);
      }

    private:

      std::vector<T> d_data;

  }; // end class

  // Define instances
  // **WARNING**  These are per particle variables.  A ParticleData object contains
  // values for all particles
  typedef MPMData<int>                     IntegerParticleData;
  typedef MPMData<double>                  DoubleParticleData;
  typedef MPMData<std::string>             StringParticleData;
  typedef MPMData<Point3D>                 Point3DParticleData;
  typedef MPMData<Vector3D>                Vector3DParticleData;
  typedef MPMData<Matrix3D>                Matrix3DParticleData;

  typedef MPMData<CellIndexVector>          VectorIntParticleData;
  typedef MPMData<CellInterpolationVector>  VectorDoubleParticleData;

  typedef MPMData<double>                  DoubleNodeData;
  typedef MPMData<Vector3D>                Vector3DNodeData;

} /* namespace BrMPM */

#endif /* MPMDATA_H_ */
