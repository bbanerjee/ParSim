/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 * MPMDataTypes.h
 *
 *  Created on: 18/10/2013
 *      Author: banerjee
 */

#ifndef MPMDATATYPES_H_
#define MPMDATATYPES_H_

#include <GeometryMath/Point3D.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>
#include <CellIndexVector.h>
#include <CellInterpolationVector.h>

#include <vector>
#include <string>
//#include <algorithm>
//#include <functional>
#include <boost/variant.hpp>

namespace BrMPM {

  // Define types
  typedef std::vector<int>          IntegerData;
  typedef std::vector<double>       DoubleData;
  typedef std::vector<std::string>  StringData;
  typedef std::vector<Point3D>      Point3DData;
  typedef std::vector<Vector3D>     Vector3DData;
  typedef std::vector<Matrix3D>     Matrix3DData;

  typedef std::vector<CellIndexVector>          VectorIntData;
  typedef std::vector<CellInterpolationVector>  VectorDoubleData;

  typedef IntegerData  IntegerParticleData;
  typedef DoubleData   DoubleParticleData;
  typedef StringData   StringParticleData;
  typedef Point3DData  Point3DParticleData;
  typedef Vector3DData Vector3DParticleData;
  typedef Matrix3DData Matrix3DParticleData;

  typedef DoubleData    DoubleNodeData;
  typedef Point3DData   Point3DNodeData;
  typedef Vector3DData  Vector3DNodeData;

  typedef VectorIntData     VectorIntParticleData;
  typedef VectorDoubleData  VectorDoubleParticleData;

  // Define the variants
  typedef boost::variant<IntegerData, DoubleData, Point3DData,
                         Vector3DData, Matrix3DData,
                         VectorIntData, VectorDoubleData> MPMVar;

  // Vistor class for zeroing out MPM data
  class ZeroVisitor : public boost::static_visitor<void>
  {
  public:
    void operator()(IntegerData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        *iter = 0;
      }
    }
    void operator()(DoubleData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        *iter = 0.0;
      }
    }
    void operator()(StringData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        *iter = "0";
      }
    }
    void operator()(Point3DData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        (*iter).set(0.0);
      }
    }
    void operator()(Vector3DData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        (*iter).set(0.0);
      }
    }
    void operator()(Matrix3DData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        (*iter).set(0.0);
      }
    }
    void operator()(VectorIntData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        (*iter).zero();
      }
    }
    void operator()(VectorDoubleData& val) const {
      for (auto iter = val.begin(); iter != val.end(); ++iter) {
        (*iter).zero();
      }
    }
  };

  // Product of two vectors
  // Assumes that operator* is available in type T
  template <typename T>
  std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    assert(vec1.size() == vec1.size());
    std::vector<T> product;
    product.reserve(vec1.size());
    for (auto iter1 = vec1.begin(), iter2 = vec2.begin() ; iter1 != vec1.end();
          ++iter1, ++iter2) {
      product.push_back((*iter1)*(*iter2));
    }
    //std::transform(vec1.begin(), vec1.end(), vec2.begin(), std::back_inserter(product),
    //               std::multiplies<T>());
    return product;
  }
}




#endif /* MPMDATATYPES_H_ */
