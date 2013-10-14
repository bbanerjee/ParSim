/*
 * MPMParticleData.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMPARTICLEDATA_H_
#define MPMPARTICLEDATA_H_

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <MPMMatrix.h>

#include <vector>
#include <string>

namespace BrMPM {

  template <class T> class MPMParticleData {

    public:

      MPMParticleData();
      virtual ~MPMParticleData();

    private:

      std::vector<T> d_data;

  }; // end class

  // Define instances
  typedef MPMParticleData<int>         IntegerParticleData;
  typedef MPMParticleData<double>      DoubleParticleData;
  typedef MPMParticleData<std::string> StringParticleData;
  typedef MPMParticleData<Point3D>     PointParticleData;
  typedef MPMParticleData<Vector3D>    VectorParticleData;
  typedef MPMParticleData<Matrix3D>    MatrixParticleData;

} /* namespace BrMPM */

#endif /* MPMPARTICLEDATA_H_ */
