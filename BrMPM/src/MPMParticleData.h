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

      MPMParticleData() {}
      virtual ~MPMParticleData() {}

      std::vector<T>::iterator begin() {return d_data.begin();}
      std::vector<T>::iterator end() {return d_data.end();}

      unsigned int size() {return d_data.size();}

    private:

      std::vector<T> d_data;

  }; // end class

  // Define instances
  // **WARNING**  These are per particle variables.  A ParticleData object contains
  // values for all particles
  typedef MPMParticleData<int>         IntegerParticleData;
  typedef MPMParticleData<double>      DoubleParticleData;
  typedef MPMParticleData<std::string> StringParticleData;
  typedef MPMParticleData<Point3D>     PointParticleData;
  typedef MPMParticleData<Vector3D>    VectorParticleData;
  typedef MPMParticleData<Matrix3D>    MatrixParticleData;
  typedef MPMParticleData<CellIndexLinear>   CellIndexLinearParticleData;
  typedef MPMParticleData<CellIndexGIMP>     CellIndexGIMPParticleData;
  typedef MPMParticleData<ShapeGradientLinear>  ShapeGradientLinearParticleData;
  typedef MPMParticleData<ShapeGradientGIMP>    ShapeGradientGIMPParticleData;

} /* namespace BrMPM */

#endif /* MPMPARTICLEDATA_H_ */
