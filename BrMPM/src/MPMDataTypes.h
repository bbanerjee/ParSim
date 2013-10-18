/*
 * MPMDataTypes.h
 *
 *  Created on: 18/10/2013
 *      Author: banerjee
 */

#ifndef MPMDATATYPES_H_
#define MPMDATATYPES_H_

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <Matrix3D.h>

#include <vector>
#include <string>
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
}




#endif /* MPMDATATYPES_H_ */
